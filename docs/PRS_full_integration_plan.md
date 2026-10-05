# Plan: Full integration of `prs_pipeline` into GDCGenomicsQC

Status: actively implemented on branch `prs_integration` — Phase 2 (single-ancestry slice + joint `multi_prscsx` + `multi_ctsleb`), Phase 3 (held-out test evaluation), and Phase 7 (every PRS rule — compute + resource-prep — runs inside the GDCGenomicsQC-owned `prs` apptainer built from `envs/prs.def`; image pushed as `oras://ghcr.io/mainsqu33ze/gdcgenomicsqc/prs:latest`) are code-complete and smoke-verified (see status + §8 session notes).
Repos involved:
- `prs_pipeline` (upstream): `/projects/standard/gdc/public/prs_methods/scripts/prs_pipeline` (git remote origin, branch `sandbox_multi_pheno` checked out)
- GDCGenomicsQC (this repo): branch `prs_integration`; active Snakefile is `workflow/Snakefile`

> **Implemented so far (phases 0–1, branch `prs_integration`):** Phase 0 baseline
> captured (fixed 5 dangling Snakefile includes; recorded pin `f71cf4f`;
> `envs/` is readable in this checkout, contrary to the earlier note). Phase 1
> (reworked 2026-09-24) **absorbs** `prs_pipeline`: the method engine is vendored
> into `workflow/scripts/prs_pipeline/` (`src/`, `templates/`, `envs/`, `LICENSE`,
> plus the transitional monolithic orchestrator) from pin `f71cf4f`; upstream is
> treated as frozen and unmaintained. `checkPRSPipelinePath` + `download_prs_resources.sh`
> now verify the vendored tree and record provenance (`prs_pipeline_path`/
> `prs_pipeline_ref`). SIF container config keys `prsMethods.containers.*` default
> to `<OUT_DIR>/containers/`. Two pre-existing blockers remain before PRS dry-runs
> go green (see Phase 0 / open questions).

## 1. Objective

Make every function the `prs_pipeline` repo offers runnable (and resumable, and
DAG-aware) as native Snakemake rules — so a single `snakemake --configfile …`
invocation can carry data prep → summary-stat generation → per-method PRS → test
evaluation. Today the workflow has partial, mostly unplugged scaffolding
(adapter rules that only validate inputs unless a `prsMethods.<method>.command`
is hand-written in the config) plus one monolithic single-ancestry rule that
shells out to `run_single_ancestry_PRS_pipeline.sh`.

## 2. What `prs_pipeline` provides (verified inventory)

### 2.1 Data prep (everything upstream of a PRS method)
| Script | Purpose |
|---|---|
| `src/genomic_preps.sh` | Attach phenotype + sex to raw PLINK data, split by ancestry (`--geno 0.05`/`--mind 0.05`), compute `.afreq` with `plink2 --freq` (required by LDpred2/lassosum2) |
| `run_prepare_prs.sh` | Orchestrates the 4-step split: `src/split_top_n_subjs.sh` → `src/run_split_plink_data.sh` → `src/generate_summary_stat_files.sh` → `src/restructure_output_dir.sh` (creates `*_gwas`, `*_study_sample`, gwas sumstats, and the split PLINK layout) |
| `run_split_data.sh` | Containerized (prsv2_latest.sif) wrapper over `src/run_split_plink_data.sh` |

### 2.2 PRS methods (all live under `src/`)
- Single-ancestry: `run_CT.sh` (C+T), `run_LDpred2.R`, `run_lassosum2.R`, `run_PRSice2.sh` (PRSice-2; requires C+T file prep), helper `prepare_sumstats.R` (bim-aligned, deduplicated 9-col sumstats), `generate_ld_matrix.R` (+ `ld_cache_dir`, `ld_matrix_dir` caching)
- Joint-ancestry: `run_PRScsx.sh` (PRS-CSx), `run_PROSPER.sh` (PROSPER, needs `singleprshelper_latest.sif`), CTSLEB (upstream CTSLEB R sources), SDPRS/SDPRX (not yet: no wrapper script in `src/`)
- Early/MVP extras: `run_viprs.sh` (VIPRS — provided-sumstats or genotype mode, needs `viprs_ref_glob` LD panels + covariate files)
- Evaluation: `score_test.sh` — scores held-out `test_sample` with training-derived parameters and writes PC-adjusted R² to `output_path/prs_pipeline/test_evaluation/`

### 2.3 Top-level orchestrators
| Script | What it does |
|---|---|
| `run_single_ancestry_PRS_pipeline.sh` | The big one. Flags `-c -l -s -P -S -B -x -v -E`; config via `-C conf`; single- or multi-phenotype mode (`summary_stats_files`+`multi_pheno_file`); runs the methods in background, then `score_test.sh` if `test_sample` set |
| `sandbox_singularity_runner.sh` | Runs the above inside `prsv2_latest.sif` with auto-bind-mount detection; submits PROSPER as a parallel SLURM job to avoid nested singularity |
| `singularity_runner.sh` | Deprecated predecessor (header says use sandbox version) |
| `run_PRS_pipeline.sh` | Older, PRS-CSx/PROSPER-centric entry; largely superseded |

### 2.4 Config & templates
- Shell-readable config files (not YAML) sourced by the scripts: `templates/single_anc_config.txt`, `templates/multi_pheno_config.txt`, `templates/1000G_simulation_config.txt`, `templates/PRSice2_sample_setup.txt`
- Containers (not in repo, pulled from GHCR): `prsv2_latest.sif`, `singleprshelper_latest.sif` — `apptainer pull oras://ghcr.io/mainsqu33ze/gdcgenomicsqc/{prsv2,singleprshelper}:latest`
- Docs: `docs/replication_workflow.md` (end-to-end ABCD/UKB replication), `docs/sandbox_singularity_runner.md`
- `envs/` exists (readable in this checkout, though `drwxrws---`, `gdc`-group-restricted) — contains `.def`/`.yml` pairs for `singlePRS`, `multiPRS`, `singlePRSHelper`; dependency list also implied by the containers

### 2.5 Key format contracts to preserve
- Single-ancestry sumstats after `prepare_sumstats.R`: `SNP CHR BP A1 A2 beta beta_se P n_eff`
- `afreq_file` (PLINK2 `.afreq`) is required for LDpred2/lassosum2
- Multi-ancestry inputs must already be *split by ancestry* and have a per-ancestry PCA/eigenvec file for C+T

## 3. Current integration in GDCGenomicsQC (verified) — what exists and what's missing

### Already wired
1. `workflow/rules/preparePRSInputs.smk` — `preparePRSInputs` rule runs `workflow/scripts/prepare_prs_inputs.sh` to build the standard split layout (`gwas/`, `anc1_plink_files/`, `anc2_plink_files/`, `metadata/`, `prs_inputs.env`, `prs_prscsx_generated.conf`, `prs_single_ancestry_<ANC1>_generated.conf`) from `phenotypeSimulation` outputs. It also writes the shared `prs_inputs.env` consumed downstream.
2. `runSingleAncestryPRS` (same file) — monolithic rule that calls the upstream `run_single_ancestry_PRS_pipeline.sh -c -l -s -P -C <generated.conf>` and touches `single_ancestry_<ANC1>.done`. **This is the only rule that actually executes a prs_pipeline method end-to-end.**
3. `workflow/rules/prsPipelines.smk` — `preparePRSMethodResources` (runs `download_prs_resources.sh` → `resources.ready`), plus 9 per-method rules (`single_ct|prsice|prscs|ldpred2|lassosum2`, `multi_ctsleb|prscsx|ldpred2|sdprs`) that all funnel through `workflow/scripts/run_prs_pipeline_adapter.sh`. The adapter validates inputs, writes `manifest.tsv`, and runs `PRS_METHOD_COMMAND` **only if the config sets one**; otherwise it writes `READY_NO_COMMAND`. → **Mostly scaffolding, unplugged by default.**
4. `workflow/rules/prosper_multi_prs_rule.smk` — real PROSPER logic (`preparePROSPERSumstats`, `preparePROSPERTuningTestingSets`, `preparePROSPERLassosumParams`, `runMultiAncestryPROSPER`) but with hard-coded `/scratch.global/saonli/...` paths and a `ref_bim` output inside the resource dir.
5. Snakefile targets: `run_preparePRSInputs`, `run_singleAncestryPRS`, `run_singleAncestryPRSPipelines`, `run_multiAncestryPRSPipelines`, `run_allPRSPipelines` (waits on all 10 `method_runs/*.done`).

### Gaps to close for "full functionality"
- [G1] No real per-method rules: 9/10 methods are stubs needing `prsMethods.<method>.command`.
- [G2] No VIPRS rule at all; no `score_test.sh` (test evaluation) rule; no multi-phenotype support in the workflow config.
- [G3] Data-prep parity: dependency on `phenotypeSimulation` outputs only. No rules mirroring `genomic_preps.sh` (real-data pheno/sex attach + ancestry split + `.afreq`) or the 4-step `run_prepare_prs.sh` splitter (for real/raw data paths like ABCD).
- [G4] Config plumbing: `prsPipeline`/`prsMethods` keys are NOT in `config/config.schema.yaml` (so `--validate` ignores them), and no example/project config defines them.
- [G5] ~~Container mismatch~~ **CLOSED 2026-09-28** (Phase 7): all PRS compute rules now declare `container: "oras://ghcr.io/mainsqu33ze/gdcgenomicsqc/prsv2:latest"` (+ `conda: ../../envs/prs.yml`); only the single `prsv2` image is used. `runSingleAncestryPRS` is still only the monolithic convenience entry — its inner `runSingleAncestryPRS` rule invokes the vendored host shell script; the per-method PRS rules are the containerized path.
- [G6] Wrapper scripts have hard-coded MSI/`saonli` paths: `run_ctsleb_wrapper.R` (`ctsleb_src`), `run_prosper_wrapper.R` (`prosper_bin`), `prosper_multi_prs_rule.smk` (multiple).
- [G7] No smoke/verification path for the PRS stages; existing `.SLURM`/configs don't exercise them.

## 4. Design decisions to confirm before coding

These change how every later phase is written. Get sign-off first.

1. **Integration style per method.** Options:
   - (A) *Native rules* wrapping each `src/*.sh|R` script with explicit inputs/outputs (recommended — matches Snakemake, resumable, per-method `*.done` already anticipated by the DAG).
   - (B) Keep the *monolithic* `run_single_ancestry_PRS_pipeline.sh` as one rule per ancestry (only what exists today).
   - (C) Hybrid: native rules for everything, keep `runSingleAncestryPRS` as a convenience that resolves to the group of method rules.
   Recommend (C). If (A) is chosen, prefer wrapping upstream scripts rather than re-implementing their internals.
2. **Vendor vs. external dependency.** Does GDCGenomicsQC (a) treat `prs_pipeline` as an external, pinned dependency located by config (`path_repo`), or (b) copy/port the needed scripts into `workflow/scripts/`? Recommend (a) + a lockstep commit pin, because `prs_pipeline` is actively developed (see branch `snakemake_prep` which already trims it for snakemake consumption).

   > **DECIDED 2026-09-24: (b) absorb/vendor.** `prs_pipeline` will **not** be
   > maintained or updated; it is being absorbed into GDCGenomicsQC. The method
   > engine is vendored at `workflow/scripts/prs_pipeline/` (source pin
   > `f71cf4f`, provenance in `VENDORED.md`), and Phase 2 replaces the monolithic
   > orchestrator with native per-method rules. `prs_pipeline_path` remains as an
   > optional override for legacy configs. Prefer rules; drop to vendored scripts
   > only where a rule is impractical.
3. **Container strategy.** Options: (i) run method rules inside `prsv2_latest.sif`/`singleprshelper_latest.sif` (what upstream uses), (ii) build GDCGenomicsQC-owned images from the same `.def`/env definitions, or (iii) per-rule conda envs for the R/Python-only methods. Recommend (i) with bind mounts for `REF`/`PRS_RESOURCE_DIR`, exposing the SIF path as config; revisit if `singleprshelper` becomes unavailable.
4. **Output layout.** Keep the current `prs_inputs/<ANC1>_<ANC2>/method_runs/` + `*.done` convention, and map upstream's `output_path/prs_pipeline/<Method>/` outputs underneath each method's `out_dir`? (Recommend yes — this is what `run_prs_pipeline_adapter.sh --out-dir` already expects.)

## 5. Phased steps

### Phase 0 — Baseline & freeze
- [x] Diff `/…/prs_pipeline` `snakemake_prep` branch vs `sandbox_multi_pheno`: `snakemake_prep` (b7c94fd) deletes ~4350 lines / 18 files — removes `score_test.sh`, `run_PROSPER.sh`, `generate_ld_matrix.R`, `genomic_preps.sh`, `run_split_data.sh`, `sandbox_singularity_runner.sh`, multi-phenotype dispatch, and PRS-CSx/VIPRS/PROSPER flags. **Decision:** base integration on `sandbox_multi_pheno` (full-featured; `snakemake_prep` is too trimmed for score-test + multi-ancestry phases).
- [x] Record the pinned commit hash of `prs_pipeline` used for integration: **`f71cf4f` (`sandbox_multi_pheno`, 2026-08-20, `f71cf4f1031ff6231bc9f8abc9f96ebdc6f17bdf`)**. Note: `envs/` **is** readable in this checkout (`drwxrws---`, readable by `gdc` group members) — the `.def`/`.yml` pairs for `singlePRS`, `multiPRS`, `singlePRSHelper` are visible, so container sources no longer need to be inferred.
- [x] Dry-run current targets to capture the present baseline. **Requires fixing 5 dangling `include:`s first**: `workflow/Snakefile` referenced `abcdBedBridge.smk`, `longitudinal_hm3_heritability.smk`, `longitudinal_hm3_spline_simulation.smk`, `longitudinal_hm3_spline_heritability.smk`, `snpHerit_abcd_adjhe_by_ancestry.smk` — never committed in any branch. Removed those 5 include lines (`rule all`/PRS rules don't reference them). After the fix, the baseline is:
  - `snakemake -n run_allPRSPipelines --configfile config/msiToy.yaml` → fails at `preparePROSPERTuningTestingSets` (`ProtectedOutputException`): it writes `ref_bim.txt` into the shared resource dir `/scratch.global/saonli/GDCGenomicsQC/prs_resources/software/PROSPER/` where a stale, non-writable copy already exists (G6 hard-coded path; Snakemake 9 raises on non-writable existing outputs). Fix belongs to Phase 2 (PROSPER rework) — the output should be a per-run file under `PRS_RESOURCE_DIR` or `PRS_METHOD_RUN_DIR`.
  - `snakemake -n run_singleAncestryPRS` / `run_singleAncestryPRSPipelines` → fails earlier at `AmbiguousRuleException` for `AFR/initialFilter.pgen` (`convertPlinkSingleFile` vs `mergeChromosomesAndFilter` can both produce it during ancestry-subset conversion). Pre-existing QC topologies/`ruleorder` issue, unrelated to PRS, blocks any end-to-end PRS dry-run until resolved (relates to Phase 9 verification).
  - Also: `--validate` no longer exists in the installed Snakemake 9.20.0 (`/projects/standard/gdc/public/envs/snakemake/bin`); `workflow/Snakefile` imports `validate` from `snakemake.utils` but never calls it, so schema validation is currently inert regardless of phase-8 schema work.

### Phase 1 — Absorb `prs_pipeline` (vendored engine + provenance)
- [x] **Vendor the method engine.** Copy `src/` (25 method/data-prep scripts incl. `PRSice.R` + `PRSice_linux`), `templates/`, `envs/` (.def/.yaml), and `LICENSE` into `workflow/scripts/prs_pipeline/` via `git archive` of pin `f71cf4f` (exact commit content, excludes the 10 GB of untracked SIFs/`temp_ld`/logs). Added `VENDORED.md` provenance. `workflow/scripts/prs_pipeline/run_single_ancestry_PRS_pipeline.sh` kept solely as the transitional fallback for `runSingleAncestryPRS` until Phase 2.
- [x] Repoint defaults at the vendored tree: `PRS_PIPELINE_PATH` in `preparePRSInputs.smk` now resolves to `workflow/scripts/prs_pipeline` (config `prs_pipeline_path` is only an optional override); `runSingleAncestryPRS` default script and `prepare_prs_inputs.sh`'s standalone default follow. `prs_pipeline_ref` (default `f71cf4f…`) is recorded **as provenance**, not as a runtime checkout check.
- [x] `checkPRSPipelinePath` (`preparePRSInputs.smk`) now verifies key vendored files exist and records `prs_pipeline_path/ref/vendored=true/checked_at` to `<OUT_DIR>/.prs_pipeline_path.checked`. No git HEAD logic (the vendored dir's git HEAD is GDCGenomicsQC's, not upstream's).
- [x] `download_prs_resources.sh`: `--prs-pipeline-dir` / `--prs-pipeline-ref` / `--prs-pipeline-sif` / `--prs-helper-sif` record provenance in `resources.ready`; file-existence checks always run; git HEAD-vs-pin verification only when the configured path is a genuine upstream checkout (origin remote contains `prs_pipeline`) — otherwise it correctly reports the vendored copy. Also hardened arg parsing (skip `--flag`-like values) so empty optional args no longer swallow subsequent flags.
- [x] Add the two SIF paths as config (`prsMethods.containers.prsv2`, `prsMethods.containers.singleprshelper`); defaults resolve to `<OUT_DIR>/containers/<name>` (SIFs are **not** vendored); GHCR pull commands documented in code + `example_config.yaml` + schema. Recorded in `resources.ready` (consumed in Phase 7).

### Phase 2 — Replace 9 stub methods with real per-method rules
Status (2026-09-24): **all five single-ancestry methods (`single_ct`, `single_ldpred2`, `single_lassosum2`, `single_prsice`, `single_prscs`) are now native rules** invoking the vendored scripts, with shared plumbing rules. Verified by smoke tests: CT/LDpred2/lassosum2 on fabricated toy data; PRSice2 on the same; PRS-CS on **real chr22 HapMap3 SNPs** (needs SNP overlap with the LD reference, so fabricated toys are unusable for it). `run_singleAncestryPRSPipelines` resolves with only the multi-ancestry stubs unsatisfied.
- [x] **Single-ancestry alignment as its own rule**: `alignSumstatsForPRS` = `prepare_sumstats.R --input <target_single> --bim <study>.bim --n_total <n> --output …/gwas/CT_PRSice2_summary_stat_file.txt` (this is exactly what the monolithic rule does internally). `n_total_gwas` configurable via `prsPipeline.n_total_gwas` (schema default 31968).
- [x] **Study-sample phenotype rule**: `makeStudyPhenoFile` (awk from `study.fam` col 6) + **study-bed→RDS rule**: `convertStudyBedToRDS` (`snp_readBed`), with stale `.rds`/`.bk` cleanup.
- [x] **LD matrix rule**: `generateLDMatrix` wraps `generate_ld_matrix.R` (map.rds, g_idx.rds, chr*_corr.rds) into `prsPipeline.ld_matrix_dir` (default `<PRS_OUT_DIR>/ld_matrix/<ANC1>`); per-method `prsMethods.<m>.ld_matrix_dir` overrides via `prs_method_ld_matrix_dir()`. Rule-conditional `input.ld_map` keeps both shared and external matrix dirs wired.
- [x] **`single_prsice`**: native rule wrapping vendored `src/run_PRSice2.sh` (PRSice.R wrapper over the vendored `PRSice_linux` binary + plink clump). `run_PRSice2.sh` accepts `PRSICE_CMD`/`PLINK_CMD` env overrides (backward-compatible); rule passes `Rscript PRSice.R --prsice <binary>` and the configured `path_plink`. `prsPipeline.binary_target` (default `"F"`) controls `--binary-target`. Verified end-to-end on the fabricated toy (`.prsice`, `.best`, `.snps` + clumped set + plots).
- [x] **`single_prscs`**: no standalone PRS-CS script existed in the engine, so a new vendored wrapper `src/run_PRScs.sh` runs **PRScsx.py in single-population mode** (= PRS-CS for that ancestry), then combines per-chr weights, `plink2 --score`, and `PRS_sscore_to_R2.R`. Config via `prsMethods.single_prscs.{path_code, ld_ref_dir, seed, path_python}` (defaults under `PRS_RESOURCE_DIR/software/PRScsx`, `…/ld/prs_csx/ref` from `download_prs_resources.sh`). Requires a python with numpy+scipy+h5py (prsv2 container or `path_python`). Verified end-to-end through Snakemake on **real chr22 HapMap3 SNPs** (from `snpinfo_mult_1kg_hm3`) fabricated into a toy study: `preparePRSInputs → align → runSingleAncestryPRSCS` produce combined weights (265 SNPs kept), `.sscore`, and R² (note: sst reformat maps aligned cols `$1,$4,$5,$6,$7`, NOT the joint script's `$1,$3,$4,$5,$6`).
- [x] **`multi_prscsx`**: native rule `runMultiAncestryPRSCSx` replacing the adapter stub — writes a KEY=VALUE temp config and calls the patched vendored `src/run_PRScsx.sh` (joint multi-pop mode of PRScsx.py: `--sst_file=<training>,<target> --n_gwas=<N_train>,<N_target> --pop=<anc2>,<anc1>` over the 8-col `target_sumstats.txt`/`training_sumstats.txt`), producing `PRScsx_<ANC>_combined_weights.txt` + `PRScsx_joint_<ANC>_score.sscore` + `{ANC1,ANC2}_PRS_sscore_Rsqr.txt` for **both** ancestries plus `multi_prscsx.done`. Config via `prsMethods.multi_prscsx.{path_code, ld_ref_dir, seed, path_python}` (defaults fall back to the `single_prscs` values). **Smoke-verified 2026-09-28** end-to-end on the real chr22 HapMap3 toys (joint 2-pop run; AFR R² 0.69 / EUR R² 0.78 toy values).
- [x] **`multi_ctsleb`**: native rule `runMultiAncestryCTSLEB` (in `workflow/rules/ctsleb_multi_prs_rule.smk`, now included in the Snakefile; stub removed from `prsPipelines.smk`) drives `workflow/scripts/run_ctsleb_wrapper.R` against the CTSLEB R sources **vendored at `workflow/scripts/CTSLEB/R`** (30 files + LICENSE + PROVENANCE; overridable via `prsMethods.multi_ctsleb.software_dir`, e.g. the `download_prs_resources.sh` `ctsleb.tar.gz`). The wrapper runs the full CT-SLEB pipeline (`SetParamsFarm`→`dimCT`→best-SNP tuning R²→`CalculateEBEffectSize`→`PRS_Clean`→SuperLearner `SL.glmnet`+`SL.ridge`→validation R² →`ExtractFinalBetas`), re-aligning plink2 `.sscore` rows to the study `.fam` order. Inputs from the existing `prepareCTSLEBSumstats` (8-col target/training sumstats + BIM→BP/rs_id) and `prepareCTSLEBPhenotypes` (study-fam col6 split into tuning/validation halves). LD-reference PLINK prefixes default to the ANC1/ANC2 study samples; point `prsMethods.multi_ctsleb.{target_ref_plink, aux_ref_plink}` at 1000G panels for real analyses. Outputs under `method_runs/multi_ctsleb/`: `{ANC1}_ctsleb_{sl_validation_r2.txt, final_coefficients.txt, dimct_prs.tsv, eb_prs.tsv, tune/validation_clean.tsv, best_snps.tsv}` + `multi_ctsleb.done` (validation is the **ANC2-as-reference → ANC1 target** prediction; toy R² 0.98). Note: `prsPipeline.generated_input_dir` is the **PRS output root** (`PRS_OUT_DIR`) — CTSLEB's sumstat/pheno prep rules key off the same key, so it must point at the prs_inputs output dir, not a subdir.
- [ ] **Multi-ancestry methods** `multi_prosper`, `multi_sdprs` (`multi_ldpred2` is also still a stub in `prsPipelines.smk`): replace adapter-only shells with native rules; PROSPER additionally needs the `ref_bim.txt` relocation listed under Blockers.
- [ ] **`afreq_file` rule**: run `plink2 --freq` on the study PLINK set to emit `.afreq` as an output consumed by LDpred2/lassosum2 (currently passed through from `prsPipeline.afreq_file` config only; schema already defines the key).
- [ ] Remove G6 hard-coded paths in `run_prosper_wrapper.R` by sourcing them from `prs_inputs.env` / config instead of literals (`run_ctsleb_wrapper.R` no longer hardcodes a path — it takes `--ctsleb-src` defaulting to the vendored `workflow/scripts/CTSLEB/R`).
- Maintenance note: vendored `src/run_lassosum2.R` was patched for compatibility with bigsnpr ≥ 1.12 (`snp_match` first-table index column is `_NUM_ID_.ss`, not `_NUM_ID_1`); `src/run_PRSice2.sh` gained `PRSICE_CMD`/`PLINK_CMD` env overrides (backward-compatible); `src/run_PRScsx.sh` gained `set -euo pipefail` + `seed`/`path_python`/`rscript` config vars (the joint-mode sst reformat `$1,$3,$4,$5,$6` is correct for our 8-col `target_sumstats.txt`); new vendored `src/run_PRScs.sh` added for single-ancestry PRS-CS (no upstream single-PRS-CS script existed). Newly vendored CTSLEB sources (`workflow/scripts/CTSLEB/R`) had two small robustness patches in `PRS_Clean.R` (from the G6 checkout): `as.data.frame()` coercions + `drop = FALSE` so single-surviving-PRS-column runs don't collapse, and `Tune_Y[[1]]` use so a 1-col pheno frame (as the tutorial passes) works. `VENDORED.md`/`PROVENANCE.md` record vendoring policy — the repo is now the maintained source of truth for the absorbed engine.
- Smoke harness: `workflow/scripts/make_smoke_toy.py` fabricated REAL chr22 HapMap3 bfiles with a phenotype correlated to 2 'causal' SNPs (guarantees real glm p<1e-3 so every C+T threshold bin is populated) and `workflow/scripts/smoke_single_ancestry.sh`, a verbose SLURM-submittable end-to-end test (`sbatch workflow/scripts/smoke_single_ancestry.sh`) that builds toys, seeds the sim subgraph, dry-runs + runs `run_singleAncestryPRSPipelines` (+ the `multi_prscsx` slice via `runMultiAncestryPRSCSx` and the `multi_ctsleb` slice via `runMultiAncestryCTSLEB`), and verifies all seven `*.done` markers + key outputs. Verified: fresh clean run in ~897 s, all 7 methods PASS (single ×5 + multi-prscsx + multi-ctsleb, plus `scoreTestPRS`).

### Phase 3 — Test evaluation (`score_test.sh`)
Status (2026-09-28): **implemented and smoke-verified** — `scoreTestPRS` bridges the native per-method `method_runs/` layout to the upstream `{CT,LDpred2,lassosum2,PRSice2,PRScsx}` layout via a `test_evaluation/_stage` symlink dir, provisions `PRSice`+plink on PATH (correcting for `score_test.sh`'s `set -eu` + bare-name calls), runs vendored `src/score_test.sh` with all `--ran-* true`, and `mv`s `<label>_{results,scores}.txt` into `method_runs/test_evaluation/`. Gated on `prsPipeline.test_sample`.
- [x] Rule `scoreTestPRS` takes `method_runs/` outputs + a configured `test_sample` PLINK prefix and runs `src/score_test.sh … --path-repo …` producing `<method>_results.txt`/`<method>_scores.txt` under `method_runs/test_evaluation/`.
- [x] Config keys `prsPipeline.test_sample`/`test_pca_eigenvec_file` (mirroring upstream); the rule builds only when set (`[rules.scoreTestPRS.output.done if … else []][0]` in `run_singleAncestryPRSPipelines`/`run_allPRSPipelines`); new target `run_scoreTestPRS`.

### Phase 4 — Multi-phenotype support
- [ ] Support `prsPipeline.summary_stats_files` (comma-separated) + `prsPipeline.multi_pheno_file` OR a workflow-native `expand()` over a phenotype list, generating per-phenotype `method_runs/<pheno>/<method>.done` (upstream writes `output_path/prs_pipeline/<pheno>/…`). Decide in Phase 0 whether multi-pheno is in-scope for the first cut.

### Phase 5 — Data-prep parity for real (non-simulated) data
- [ ] Add optional rules mirroring `genomic_preps.sh`: attach phenotype/sex, split by ancestry keep-list, `plink2 --freq` → `.afreq`, feeding the same `preparePRSInputs` inputs the simulation path currently produces.
- [ ] Optionally add the 4-step splitter (`split_top_n_subjs`→`run_split_plink_data`→`generate_summary_stat_files`→`restructure_output_dir`) only if real-data runs need it; otherwise document that `run_prepare_prs.sh` remains the external pre-step for that path.
- [ ] Generalize `prepare_prs_inputs.sh` so the PLINK sources come from config (real-data dirs or per-ancestry `unrelated.pgen` from QC) rather than only `phenotypeSimulation` outputs.

### Phase 6 — VIPRS
- [ ] Add `rule runSingleAncestryVIPRS` wrapping `src/run_viprs.sh --c <temp config>` with `provided_sumstats` (recommended; no GWAS genotypes needed) or `bfile_gwas_input`+`covariate_file_gwas` modes, `covariate_file_study_sample`, and `viprs_ref_glob` LD panels.
- [ ] Add `viprs_ref_glob` handling to `download_prs_resources.sh` (VIPRS LD panels are third-party; document download per VIPRS docs) and config keys under `prsMethods.single_viprs`.

### Phase 7 — Containers & environments
- [x] **Run every PRS compute rule inside the upstream `prsv2` container** (`oras://ghcr.io/mainsqu33ze/gdcgenomicsqc/prsv2:latest`, per Phase 0 decision (i)). Added `container:` + `conda: "../../envs/prs.yml"` (fresh `envs/prs.yml` parity recipe; no CTSLEB GitHub install — the wrapper sources the vendored `workflow/scripts/CTSLEB/R`) to all 13 compute rules: `preparePRSInputs`, `alignSumstatsForPRS`, `convertStudyBedToRDS`, `generateLDMatrix`, `runSingleAncestry{CT,LDpred2,Lassosum2,PRSice,PRSCS}`, `runSingleAncestryPRS` (monolithic), `scoreTestPRS`, `runMultiAncestryPRSCSx`, `runMultiAncestryCTSLEB`. Apptainer binds (`/scratch.global,/projects,/spaces,/home`) match the production profiles; images cache in the gitignored `workflow/.snakemake/singularity/`.
- [x] **Interpreter resolution**: rules no longer resolve `shutil.which("python3")/Rscript` on the host (wrong paths inside a container) — `path_python` defaults to the bare `python` on the image PATH (`PRS-CS`/`PRS-CSx` endpoints fall back to `python` when unset) and `rscript` is the bare image `Rscript`. Host `plink`/`plink2`/python config keys are no longer set by the smoke.
- [x] **`runSingleAncestryPRS` mono-containerized**; its PRSice step drives the vendored `PRSice.R`+`PRSice_linux` via `PRSICE_CMD` (the published `prsv2` image has no `PRSice` command).
- [x] **`envs/prs.def` now bakes a self-contained image**: copies the vendored `prs_pipeline` (incl. static `PRSice_linux`) and `CTSLEB` into `/opt/gdcgenomicsqc/`; clones the upstream method repos at pinned commits into `/opt/gdcgenomicsqc/software/` (PRScs `5330390`, **PRScsx `9a77968` — the smoke-tested clone HEAD**, CTSLEB `4e37846`, PROSPER `b4b4a1b`, SDPRX `1fd89bc`, all frozen 2026-09-28, `.git` stripped); installs `git`/`curl`/`ca-certificates` so the resource-prep step *could* also run in-container.
- [x] Keep the following **host-side** (documented exceptions): `preparePRSMethodResources`/`checkPRSPipelinePath` (resource setup/file checks; git/curl now baked into the image, so these may move in-container once a GDC-owned build is pushed — the published `prsv2` image still lacks git/curl), `makeStudyPhenoFile`/`prepareCTSLEB{Phenotypes,Sumstats}` (pure awk/`run:`-block), `runMultiAncestry{LDpred2,SDPRS}` (stubs).
- [ ] Nested-container trap still open: PROSPER inside a container → parallel SLURM job (as upstream `sandbox_singularity_runner.sh`); must be handled or explicitly rejected when `multi_prosper` is wired.
- [ ] Build + push a GDCGenomicsQC-owned image from `envs/prs.def` (`cd envs && apptainer build --fakeroot prs.sif prs.def` → push to GHCR), then point the rules at it and flip `preparePRSMethodResources` in-container.

### Phase 8 — Config schema, docs, and entry targets
- [ ] Add `prsPipeline` and `prsMethods` blocks to `config/config.schema.yaml` (all keys currently read via `config.get(...)` in `prsPipelines.smk`, `preparePRSInputs.smk`, `prosper_multi_prs_rule.smk`: `generated_input_dir`, `single_ancestry_script/flags`, `phenotype_index`, `gwas_fraction`, `seed`, `path_plink2`, `resource_dir`, `download_software`, per-method `command/ld_ref_dir/ld_ref_prefix/ld_matrix_dir/software_dir`, `plink2`, `containers`).
- [ ] Add a `prsMethods`+`prsPipeline` example block to `config/example_config.yaml` and to one toy config (`msiToy.yaml`/`sandboxToy.yaml`) so `run_allPRSPipelines` is smoke-testable.
- [ ] Update Snakefile `rule all` / stage convenience targets to include `run_allPRSPipelines`-style inputs only when PRS config is present (guarded), and document the `run_*` targets in README/`docs/usage.rst`.
- [ ] Update `AGENTS.md` PRS section and add a short `docs/prs_methods.md` mapping upstream scripts → workflow rules.

### Phase 9 — Verification & hardening
- [ ] Config validation: `snakemake --validate --configfile config/msiToy.yaml` passes with PRS keys.
- [ ] Dry runs: `snakemake -n run_singleAncestryPRS`, `-n run_allPRSPipelines` resolve without unmet input errors.
- [ ] Smoke run on toy data (AFR+EUR, small `phenotypeSimulation`): end-to-end from `preparePRSInputs` through all `.done` markers; confirm sumstats column contract (`SNP A1 A2 BETA SE P N`), `.afreq` presence, and `score_test` R² files.
- [ ] Confirm each method rule is idempotent/resumable via its `*.done` marker and logs under `OUT_DIR/logs/`.
- [ ] Run `git status`/diff review to ensure no stale backups (`#…#`, `~`, `.bak`, `.before_*`) got edited in `workflow/rules/`.

## 6. Suggested first implementation slice (smallest end-to-end)

Do not attempt all phases at once. Suggested order for a first working cut:
1. Phase 1 (absorb/vendor the engine + provenance checks).
2. ~~Phase 2 for `single_ct` + `single_ldpred2` + `single_lassosum2` only (covers the shared `prepare_sumstats` + `afreq` + LD-matrix plumbing and produces real `*.done`)~~ **DONE 2026-09-24** (smoke-tested end-to-end on fabricated toy data; see Phase 2 status above).
3. ~~Phase 2 remainder of the single-ancestry slice: `single_prsice` (vendored `run_PRSice2.sh` + `PRSice.R`/`PRSice_linux`) and `single_prscs` (new `src/run_PRScs.sh`, PRScsx.py single-pop)~~ **DONE 2026-09-24** (PRSice smoke on fabricated toy; PRS-CS smoke on a real chr22 HapMap3 toy through Snakemake).
4. ~~Phase 3 (`score_test.sh`) so single-ancestry is evaluable~~ **DONE 2026-09-28** (evaluable: `scoreTestPRS` produces `<method>_{results,scores}.txt` per method on a held-out sample).
5. The remaining multi-ancestry methods (Phase 2 remainder: `multi_prosper`, `multi_sdprs`) — plus wire `multi_ldpred2` (still a stub in `prsPipelines.smk`) — followed by VIPRS (Phase 6), multi-pheno (Phase 4), and data-prep parity (Phase 5) as separate increments.
6. ~~Phase 7 (containers)~~ **DONE 2026-09-28** (every PRS compute rule runs inside `oras://ghcr.io/mainsqu33ze/gdcgenomicsqc/prsv2:latest` via `--use-apptainer`, smoke-verified; `envs/prs.yml`+`prs.def` added).
7. Phase 8 (schema/example config) throughout, since each increment adds config keys.

## 7. Open questions for the maintainer (blocking Phase 0 sign-off)

Resolved during Phase 0/1:
- Integration style: native per-method rules (Phase 2), keeping `runSingleAncestryPRS` as the group convenience — hybrid (C). Prefer rules; drop to vendored scripts only where a rule is impractical.
- Absorption: `prs_pipeline` is vendored into `workflow/scripts/prs_pipeline/` at pin **`f71cf4f` (`sandbox_multi_pheno`)**; upstream is frozen/unmaintained. `prs_pipeline_path`/`prs_pipeline_ref` config keys record provenance (Phase 1, reworked 2026-09-24).
- Container approach: run method rules inside upstream `prsv2_latest.sif`/`singleprshelper_latest.sif` (recommendation (i)); paths configurable via `prsMethods.containers.*` (Phase 1) and consumed in Phase 7. → **resolved 2026-09-28**: every PRS compute rule declares `container: "oras://ghcr.io/mainsqu33ze/gdcgenomicsqc/prsv2:latest"` + `conda: "../../envs/prs.yml"`; only the single `prsv2` image is used (the `singleprshelper` SIF was not needed — PROSPER, the only consumer, is deferred).

Still open (blocking green PRS dry-runs):
- PROSPER `preparePROSPERTuningTestingSets` writes `ref_bim.txt` into the shared `prs_resources/software/PROSPER/` dir where a stale, non-writable file exists → Phase 2 must relocate this output under `PRS_RESOURCE_DIR`/`PRS_METHOD_RUN_DIR` (and de-hardcode the `/scratch.global/saonli/...` paths).
- Pre-existing `AmbiguousRuleException` at `<ANC>/initialFilter.pgen` (`convertPlinkSingleFile` vs `mergeChromosomesAndFilter`) blocks the full upstream QC → PRS subgraph from dry-running; the Phase 2 smoke test worked around it via a pre-seeded `phenotypeSimulation.input_prefixes` toy OUT_DIR. Resolve topology/`ruleorder` so the Phase 9 smoke run can proceed on real QC output.
- Multi-phenotype mode remains out of scope for the first cut (deferred to Phase 4).
- Confirm a writable dev `OUT_DIR` + toy PRS resources for the smoke test (Phase 9). (Smoke harness for the Phase 2 slice: `phenotypeSimulation.input_prefixes` pointing at existing bfiles so the QC/simulation subgraph is skipped; test config `/tmp/opencode/smoke3.yaml`.)

## 8. Session notes

### Session 2026-09-28 (5) — all PRS compute rules containerized incl. monolithic; image recipe now bakes the software
- Per maintainer directive ("use singularity images for all dependencies; new images can be built if needed") completed the last un-containerized PRS rule: **`runSingleAncestryPRS`** (monolithic `run_single_ancestry_PRS_pipeline.sh`) now declares `container:`/`conda:` like the other 12. Since the published `prsv2` image has no `PRSice` command, its PRSice sub-method drives the vendored `PRSice.R`+`PRSice_linux` via `export PRSICE_CMD=...` (same pattern as `runSingleAncestryPRSice`). Dry-run green.
- **`envs/prs.def` downgraded to a full self-contained image recipe** (built from `envs/`):
  - `%files` bakes the vendored `workflow/scripts/prs_pipeline` (incl. static `PRSice_linux`) and `workflow/scripts/CTSLEB` into `/opt/gdcgenomicsqc/` — the image no longer needs the repo mount for the vendored software (the workflow still uses repo paths by default).
  - `%post` installs `git`/`curl`/`ca-certificates` (apt) and clones the upstream method repos at frozen commits into `/opt/gdcgenomicsqc/software/`, `.git` stripped for size: **PRScsx `9a77968755e576d58aa43ed205b9d0f02eac6468` (the exact clone HEAD the smoke-tested PRS-CS(x) runs use — not the upstream HEAD `37eaeb3`)**, PRScs `5330390…`, CTSLEB `4e37846…`, PROSPER `b4b4a1b…`, SDPRX `1fd89bc…`.
- Audit (13/13 compute rules have `container:`+`conda:`): remaining host rules are `preparePRSMethodResources`/`checkPRSPipelinePath` (setup + file checks), `makeStudyPhenoFile`/`prepareCTSLEB{Phenotypes,Sumstats}` (awk/`run:`-block), `runMultiAncestry{LDpred2,SDPRS}` (stubs). `preparePRSMethodResources` could move in-container once a GDC-owned image (with git/curl) is built+pushed — the published `prsv2` still lacks them.
- Next: build+push the GDC-owned image; then optionally point `preparePRSMethodResources` (and rule defaults like `single_prscs.path_code`) at the `/opt/gdcgenomicsqc/…` in-image paths.
- Git: **no commits made** — edits to `envs/prs.def`, `workflow/rules/preparePRSInputs.smk`, `docs/PRS_full_integration_plan.md`.

### Session 2026-09-28 (4) — PRS rules containerized (apptainer, upstream `prsv2` image)
- Per maintainer directive ("the rest of the pipeline uses singularity images … as few images as possible"), the PRS rules now run inside the upstream container `oras://ghcr.io/mainsqu33ze/gdcgenomicsqc/prsv2:latest` instead of host tooling, mirroring `workflow/rules/Relatedness.smk`'s `container:`+`conda:` pattern. Inspected the image (digest `sha256:c33e746f…`): single conda env `multiPRS` on PATH, R 4.3.0 with the required packages (bigsnpr 1.12.18, CTSLEB, SuperLearner, caret, glmnet, pROC, …), bioconda `plink`/`plink2`, conda python (numpy 1.26.4, scipy 1.17.1, h5py 3.9.0); has **no** `prsice` binary (vendored `PRSice_linux` is used), **no `git`/`wget`**.
- Added `container:` + `conda:"../../envs/prs.yml"` to all 12 compute rules: `preparePRSInputs` (its `conda:` switched off `phenotypeSim.yml`), `alignSumstatsForPRS`, `convertStudyBedToRDS`, `generateLDMatrix`, `runSingleAncestry{CT,LDpred2,Lassosum2,PRSice,PRSCS}`, `scoreTestPRS`, `runMultiAncestryPRSCSx`, `runMultiAncestryCTSLEB`. Added parity recipe `envs/prs.yml` (+ `envs/prs.def` micromamba build documenting rebuild/push; no CTSLEB GitHub install — the wrapper sources the vendored `workflow/scripts/CTSLEB/R`).
- Interpreter fixes (container-correct): PRS-CS/PRS-CSx `path_python` no longer falls back to `shutil.which("python3")` (a host path, meaningless inside a container) → empty string, so the vendored `run_PRScs*.sh` use bare `python` on the image PATH (their existing `PYTHON_CMD` fallback); `rscript` → bare `"Rscript"`. Host `path_plink`/`path_plink2` dropped from the smoke config → defaults `plink`/`plink2` resolve to the image binaries.
- **Smoke-verified 2026-09-28**: fresh full run (`GDCQC_SMOKE_WORK=/scratch.global/baron063/testing/GDCQC_apptainer_full`) PASSED in **626 s** entirely via `apptainer exec … workflow/.snakemake/singularity/<hash>.simg` with binds `/scratch.global,/projects,/spaces,/home`. All 5 singles + `multi_prscsx` + `multi_ctsleb` (AFR validation R² 0.900 on the toy) + `scoreTestPRS` green.
- Host-side exceptions (by design): `preparePRSMethodResources`/`checkPRSPipelinePath` (need `git`, absent from the image), `makeStudyPhenoFile`/`prepareCTSLEBPhenotypes`/`prepareCTSLEBSumstats` (pure awk / `run:`-block — containers don't apply to `run:`), `runMultiAncestry{LDpred2,SDPRS}` (stubs).
- Git: **no commits made** — new `envs/prs.yml` + `envs/prs.def` (untracked); edits to `workflow/rules/{prsPipelines,ctsleb_multi_prs_rule,preparePRSInputs}.smk`, `workflow/scripts/smoke_single_ancestry.sh`, `docs/PRS_full_integration_plan.md`, plus the still-uncommitted multi_ctsleb/prscsx work from sessions (2)/(3).

### Session 2026-09-28 (3) — `multi_ctsleb` (CT-SLEB; multi-ancestry reference → single-ancestry target) implemented & smoke-verified
- Vendored the CTSLEB R sources from the G6 checkout `/scratch.global/saonli/GDCGenomicsQC/prs_resources/software/CTSLEB_github/R` into `workflow/scripts/CTSLEB/R` (all 30 files byte-identical) + `LICENSE` (MIT, 2022 Haoyu Zhang) + `PROVENANCE.md`. Sources were already on the saonli scratch, so no web download was needed.
- Rewrote `workflow/scripts/run_ctsleb_wrapper.R` as a self-contained driver (previously a hardcoded-path adapter): `--ctsleb-src` resolves vendored path by default (arg → `CTSLEB_R_DIR` env → `workflow/scripts/CTSLEB/R`); runs `SetParamsFarm` → `dimCT(sum_target, sum_ref, ref_plink=ANC2-study, target_plink=ANC1-study, test_target_plink=ANC1-study)` → per-PRS tuning R² → best-SNP set → `CalculateEBEffectSize` → re-align PRS rows to study `.fam` IID order (plink2 `.sscore` order ≠ fam order) → `PRS_Clean` (SuperLearner `SL.glmnet`+`SL.ridge`, nearZeroVar filters) → held-out validation R² → `ExtractFinalBetas`. Writes `_{sl_validation_r2.txt, final_coefficients.txt, dimct_prs.tsv, eb_prs.tsv, tune_clean.tsv, validation_clean.tsv, best_snps.tsv}` + `.done`.
- New `workflow/rules/ctsleb_multi_prs_rule.smk` (native rule, now `include`d in the Snakefile; old adapter stub removed from `prsPipelines.smk`). Plumbed into the existing `prepareCTSLEBSumstats`/`prepareCTSLEBPhenotypes` rules (sumstats→CTSLEB format via BIM for BP; study-fam col6 split into `{ANC1}_HM3_tuning.pheno`/`_validation.pheno`). Config `prsMethods.multi_ctsleb.{software_dir, target_ref_plink, aux_ref_plink}` — LD-ref prefixes default to the ANC1/ANC2 study samples (the toy has no 1000G panel to test against), so real analyses must point these at 1000G.
- **Smoke-verified 2026-09-28**: full clean run (fresh workspace, `GDCQC_SMOKE_WORK=/scratch.global/baron063/testing/GDCQC_multiCTSLEB_dry`) passed in 897 s — all 5 singles + `multi_prscsx` + `multi_ctsleb` + `scoreTestPRS` green. CT-SLEB super-learner validation R² (AFR target) = 0.984 on the toy.
- Gotcha fixed mid-session: `prsPipeline.generated_input_dir` is the **PRS output root** (`PRS_OUT_DIR`), not a CTSLEB-specific input dir — the smoke config first set it to a `ctsleb_inputs` subdir, silently relocating **all** method outputs (including the previously-passing singles/prscsx) under it and failing every harness check despite successful runs. Point it at `OUT/prs_inputs/{ANC1}_{ANC2}` (its own default). Second teething bug: `PRS_Clean` died on `colSums`/`colMeans` of a 1-D vector — fixed by patching the vendored `PRS_Clean.R` and passing a 1-col data.frame to `Tune_Y`.
- Git: per maintainer policy **no commits made** — changes are working-tree (snakefile/rules, wrapper, vendored CTSLEB dir, schema/example config, smoke harness, plan doc).

### Session 2026-09-28 (2) — `multi_prscsx` (joint multi-ancestry PRS-CSx) implemented & smoke-verified
- Replaced the adapter stub `rule runMultiAncestryPRSCSx` (prsPipelines.smk) with a native rule mirroring `runSingleAncestryPRSCS`: temp `PRScsx_temp_config.txt` + patched vendored `src/run_PRScsx.sh` (added `set -euo pipefail`, `seed`/`path_python`/`rscript` config; sst reformat left as `$1,$3,$4,$5,$6`, correct for our 8-col `target_sumstats.txt`). Outputs per ancestry: `PRScsx_<ANC>_combined_weights.txt`, `PRScsx_joint_<ANC>_score.sscore`, `<ANC>_PRS_sscore_Rsqr.txt`, `<ANC>_adj_PRS_sscore_Rsqr.txt` + `multi_prscsx.done`.
- Globals `PRS_PRSCSX_{PATH_CODE,REF_DIR,SEED,PATH_PYTHON}` default from `prsMethods.multi_prscsx` → `single_prscs` → resource dir; schema + example config gained the `multi_prscsx` block (mirrors `single_prscs`).
- Verified joint call: `PRScsx.py --sst_file=<training>,<target> --n_gwas=50,50 --pop=EUR,AFR --out_name=PRScsx --seed=42`; `*** 2 discovery populations detected ***`; LD aligned across all 22 chr (only chr22 carries toy SNPs); both ANC scored; toy R² AFR 0.69 / EUR 0.78 (high because toy train/target share causal variants).
- Smoke harness now runs the multi slice too (fresh-workspace PASS ~490 s); `run_multiAncestryPRSPipelines`/`run_allPRSPipelines` still won't resolve until the other 3 multi stubs are replaced.
- Git: per maintainer policy **no commits made** — changes are working-tree (schema/example, prsPipelines.smk, vendored run_PRScsx.sh, smoke harness, plan doc).

### Session 2026-09-28 — Phase 3 (`scoreTestPRS`) implemented & smoke-verified
- Used as the briefing for next session; the 09-24 notes and §5/§6 below remain the reference for the method engine and gotchas.
- `scoreTestPRS` + `run_scoreTestPRS` + gating in `workflow/Snakefile`; config keys `prsPipeline.test_sample`/`test_pca_eigenvec_file` in schema + example config.
- Smoke harness extended: fabricates a held-out test toy (seed+200, distinct IIDs) → runs `run_scoreTestPRS` → verifies `test_evaluation/{CT,LDpred2_inf,LDpred2_grid,lassosum2,PRSice2,PRScsx_<ANC>}_{results,scores}.txt` + `scoreTestPRS.done`. Full e2e PASS in ~216 s (14/14 jobs).
- Held-out R² ≈ 0 is expected (each toy study is an independent phenotype draw) — Phase 3 verifies plumbing (weights → score → R² files), not effect size.
- Git: Phase 2 + vendored engine + harness committed as `7d5c841`; this session's Phase 3 changes are the 5 files in the working tree (config schema/example, Snakefile, prsPipelines.smk, smoke harness). `temp_test_script.sh` + `workflow/scripts/__pycache__/` are gitignored.

### Session 2026-09-24 (Phase 2)
- Branch `prs_integration`; **nothing committed** — all of Phase 2 + the smoke harness are working-tree changes. Review with `git status`/`git diff` before continuing.
- Modified: `.gitignore`, `config/config.schema.yaml`, `config/example_config.yaml`, `workflow/Snakefile`, `workflow/rules/preparePRSInputs.smk`, `workflow/rules/prsPipelines.smk`, `workflow/scripts/download_prs_resources.sh`, `workflow/scripts/prepare_prs_inputs.sh`, `docs/PRS_full_integration_plan.md`.
- New (untracked): `workflow/scripts/prs_pipeline/` (vendored engine), `workflow/scripts/make_smoke_toy.py`, `workflow/scripts/smoke_single_ancestry.sh`, `docs/PRS_full_integration_plan.md`.
- All five single-ancestry methods done + e2e smoke verified (13/13 jobs): CT, PRSice2 (best p=5e-8, R² 0.685), PRS-CS (Rsqr 0.633), LDpred2 (Grid R² 0.80), lassosum2.

### Smoke harness (how to re-run)
- `sbatch workflow/scripts/smoke_single_ancestry.sh` (or run interactively). Requires **apptainer** on PATH — every PRS rule runs in the `prsv2` container; snakemake is invoked with `--use-apptainer --apptainer-args "--bind /scratch.global,/projects,/spaces,/home"`. Images cache in the gitignored `workflow/.snakemake/singularity/` (snakemake 9.20 ignores `SNAKEMAKE_APPTAINER_PREFIX`; there is no `--apptainer-prefix` override in the harness).
- Default workspace: `/scratch.global/baron063/testing/GDCQC_PRS_integration/` (override with `GDCQC_SMOKE_WORK`). Env knobs: `GDCQC_SMOKE_KEEP=1` (reuse), `GDCQC_SMOKE_DRYRUN=1`, `GDCQC_SMOKE_CPUS`, `GDCQC_SMOKE_SNPS/SAMPLES/SEED`, `GDCQC_SMOKE_TRACE`.
- As of Session (4), the harness no longer sets `prsPipeline.path_plink/path_plink2` or `single_prscs.path_python` — inside the container they resolve from the image PATH. Full containerized run ≈ **626 s** (was ~897 s host-mode).
- Toy bfiles must (a) use REAL chr22 HapMap3 SNP IDs (PRS-CS needs overlap with `snpinfo_mult_1kg_hm3`/the LD ref) and (b) have a fam phenotype correlated with 2 "causal" SNPs, so the pipeline's real `plink2 --glm` yields p<1e-3. Without that, the smallest `run_CT.sh` q-score-range bin (0.001) is empty → `temp.0.001.profile` never written → the C+T R step dies reading the missing file.
- Ephemeral scratch/test dirs under `/tmp/opencode/` (smoke3, vendor_test, vendor_prscs, prscs_smoke, prscs_snake, smoke_dry*, smoke_e2e*) can be discarded.

### Gotchas to remember next session (each hit and solved here)
- **snakemake 9.20.0 CLI**: `--configfile` greedily consumes the following non-option tokens as *additional* config files. Put the target rule BEFORE `--configfile`, e.g. `snakemake -j 8 run_singleAncestryPRSPipelines --configfile cfg.yaml`. The form `-j 8 --configfile cfg.yaml TARGET` silently treats `TARGET` as a config path → `FileNotFoundError: 'run_singleAncestryPRSPipelines'`.
- **PATH pollution** *(now moot for the PRS rules — Session (4): they execute inside the `prsv2` container with the image's conda PATH; host R/python/plink are no longer used)*: never prepend `/projects/standard/gdc/public/envs/plink/bin` or `gdcPipeline/bin` to PATH — both ship a conda `Rscript`/`R` that lacks argparse/bigsnpr and silently shadows the working R for every rule shell. Resolve plink/plink2 to absolute paths (config `prsPipeline.path_plink/path_plink2`) instead.
- Working `Rscript` = `/common/software/install/manual/R/4.4.0-openblas-rocky8-fix/bin/Rscript` (bigsnpr 1.12.21, argparse, R.cache, …). Validate: `Rscript -e 'library(argparse); library(bigsnpr)'`. *(Still the host reference; the PRS rules use the image's R 4.3.0.)*
- Working PRS-CS python = `/projects/standard/gdc/public/envs/gdcPipeline/bin/python` (numpy 1.26.4, scipy 1.16.3, h5py 3.15.1). Rules resolve it via `prsMethods.single_prscs.path_python` (fallback `shutil.which`), never via PATH. *(In-container default: bare `python` from the image.)*
- Rules re-resolve `Rscript`/interpreter paths at DAG time (`shutil.which`) and via bare `Rscript` in shells — any PATH fix must be in effect before `snakemake` starts. *(Now mostly moot: only the non-containerized rules still do this.)*
- **Apptainer gotchas (Session (4))**: jinja `container:` URIs are resolved by snakemake against its run-dir cache `workflow/.snakemake/singularity/` (gitignored); a first `--use-apptainer` run pulls ~1.5 GB from GHCR. `apptainer exec` gets `--home <workflow-dir>` (writable under `/projects`) plus the configured binds; `preparePRSMethodResources` stays on the host because the image lacks `git`.
- `simulateBivariatePhenotypes` bfile branch hardcodes `/scratch.global/saonli/GDCGenomicsQC/CTSLEB/…` paths → the smoke bypasses the sim/QC subgraph by pre-seeding `OUT_DIR/simulations/{ANC1}_{ANC2}/*_simulation.{bed,bim,fam}` (copies of the toy bfiles, copied AFTER the toy inputs so their mtimes are newer).
- Snakemake 9.20 has **no `--validate`**; config correctness is checked via `-n` dry-runs.
- `.done` markers are intentionally 0-byte → assert with `-e`, not `-s`.

### Next increments (see §5/§6)
1. ~~Phase 3 `scoreTestPRS` (`src/score_test.sh`) — makes the single-ancestry methods evaluable~~ **DONE 2026-09-28** (smoke-verified; see §8).
2. ~~Phase 2 remainder + multi slice: `multi_prscsx`~~ **DONE 2026-09-28** (smoke-verified; see §8) and ~~`multi_ctsleb`~~ **DONE 2026-09-28** (smoke-verified; see §8). Then PROSPER (needs the `ref_bim.txt` relocation + G6 de-hardcoding), `multi_sdprs`, `multi_ldpred2`.
3. ~~Phase 7 containers~~ **DONE 2026-09-28** (all PRS compute rules in `prsv2`; smoke-verified; see §8 session (4)).
4. `afreq_file` rule (`plink2 --freq` on the study set), then G6 path cleanup in the C+T/PROSPER wrappers.

### External (non-vendored) resources the PRS-CS work references
- `/projects/standard/gdc/public/prs_methods/scripts/PRScsx` (PRScsx.py, parse_genet.py)
- `/projects/standard/gdc/public/prs_methods/ref/ref_PRScsx/1kg_ref` (`ldblk_1kg_{afr,amr,eas,eur,sas}`, `snpinfo_mult_1kg_hm3`)
- Both are wired through `prsMethods.single_prscs.{path_code, ld_ref_dir}` in the smoke config / schema defaults.