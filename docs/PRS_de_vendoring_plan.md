# Plan: De-vendor `prs_pipeline` into flat workflow scripts and image-provided tools

Status: **Implemented 2026-11-06** (sections 5.1–5.6 complete; validation in
progress). PRS full-scale fixes applied and toy smoke passed 2026-11-03.

Date: 2026-11-02

Updated: 2026-11-03 — applied section 6.6 fixes; wired the external LD matrix
(`ref/ldpred2`) into the full/reduced configs and confirmed the dry run no longer
builds `generateLDMatrix`; toy smoke passed (section 6.6.7).

Updated: 2026-11-06 — de-vendoring executed: the 12 active scripts were moved
flat into `workflow/scripts/` (carrying the §6.6.3 fixes), wrappers made
sibling-relative via `SCRIPT_DIR`, `PRS_PIPELINE_PATH`/`path_repo`/vendored
PRSice plumbing and the monolithic rule removed, `preparePRSMethodResources`
replaced by `checkPRSImageTools` (`OUT_DIR/.prs_image_tools.checked`), image tools
made authoritative in `config/ctsleb_sim_prs*.yaml`, `envs/prs.def` updated (no
vendored-tree `%files` line), and the vendored `workflow/scripts/prs_pipeline/`
tree deleted. **Deferred per maintainer:** image rebuild/push; PROSPER keeps
host `software_dir`/`plink2` overrides because `singleprshelper` lacks plink2 and
`/opt/gdcgenomicsqc`. CT-SLEB also keeps the repo `workflow/scripts/CTSLEB/R`
because the image copy predates the plink2 `.sscore` fix. Validation target:
`config/ctsleb_sim_prs_small.yaml` (§6.6.4).

## 1. Objective

Remove the runtime dependency on `workflow/scripts/prs_pipeline/` while keeping the
current PRS behavior unchanged.

Target end state:

- The 12 active PRS helper scripts live **directly** under `workflow/scripts/`
  (flat, not nested under `prs_pipeline/` or `prs_methods/`).
- PRS method tools and software come from the singularity images, not from host
  checkouts or vendored copies.
- Snakemake rules reference tools at stable image paths or via `PATH`.
- No `PRS_PIPELINE_PATH`, `prs_pipeline_path`, `prs_pipeline_ref`,
  `path_repo`, or vendored `PRSice.R`/`PRSice_linux` remains.
- Existing outputs and method behavior are preserved.

## 2. Scope and assumptions

### In scope

- Move the active vendored scripts into `workflow/scripts/`.
- Rewire all active Snakemake rules and shell/R wrappers.
- Delete the vendored `prs_pipeline` tree and its provenance/config plumbing.
- Use image-provided PRSice, plink/plink2, PRScsx, PRScs, CT-SLEB, and PROSPER.
- Remove host-path overrides from example/smoke configs where they duplicate image
  content.
- Update the PRS image definition and verification/smoke tooling.
- Validate with dry runs, the smoke workflow, and the full `ctsleb_sim_prs.yaml`
  workflow.

### Out of scope unless requested

- Reimplementing PRS method math.
- Completing `multi_ldpred2`, `multi_sdprs`, VIPRS, multi-phenotype, or real-data
  split parity.
- Baking multi-GB LD reference panels into the SIF.

### Open decision

LD reference panels are assumed to remain read-only data supplied through
`ld_ref_dir` and bind-mounted from `/projects` or `/scratch`. If they must also be
baked into the singularity image, add 1–3 days and account for a much larger image.

## 3. Scripts to move

Move these from `workflow/scripts/prs_pipeline/src/` to `workflow/scripts/`:

| Current path | New path |
|---|---|
| `src/prepare_sumstats.R` | `workflow/scripts/prepare_sumstats.R` |
| `src/generate_ld_matrix.R` | `workflow/scripts/generate_ld_matrix.R` |
| `src/run_CT.sh` | `workflow/scripts/run_CT.sh` |
| `src/CT_script.R` | `workflow/scripts/CT_script.R` |
| `src/run_LDpred2.R` | `workflow/scripts/run_LDpred2.R` |
| `src/run_lassosum2.R` | `workflow/scripts/run_lassosum2.R` |
| `src/run_PRSice2.sh` | `workflow/scripts/run_PRSice2.sh` |
| `src/run_PRScs.sh` | `workflow/scripts/run_PRScs.sh` |
| `src/run_PRScsx.sh` | `workflow/scripts/run_PRScsx.sh` |
| `src/PRS_sscore_to_R2.R` | `workflow/scripts/PRS_sscore_to_R2.R` |
| `src/score_test.sh` | `workflow/scripts/score_test.sh` |
| `src/evaluate_test.R` | `workflow/scripts/evaluate_test.R` |

No filename collisions exist with the current `workflow/scripts/` contents.
Keep the existing filenames to minimize diffs.

Do **not** move `PRSice.R` or `PRSice_linux`. Both are third-party; use the
image-provided `PRSice` command instead.

## 4. Current active dependency map

- `workflow/rules/prsPipelines.smk:36` defines
  `PRS_SRC = PRS_PIPELINE_PATH / "src"`.
- `workflow/rules/preparePRSInputs.smk` defines and validates
  `PRS_PIPELINE_PATH` at `:7-12` and `:51-83`.
- `workflow/rules/preparePRSInputs.smk:136` passes `prs_pipeline_dir` into
  `prepare_prs_inputs.sh`.
- `workflow/rules/preparePRSInputs.smk:180` references the monolithic script.
- `workflow/rules/preparePRSInputs.smk:200` builds `PRSICE_CMD` from vendored
  `PRSice.R` + `PRSice_linux`.
- `workflow/rules/prsPipelines.smk:318`, `:554`, `:832` write
  `path_prs_pipeline`/`prs_pipeline` into generated configs.
- `workflow/rules/prsPipelines.smk:479-488` uses vendored PRSice.
- `workflow/rules/prsPipelines.smk:619-652` uses vendored PRSice and
  `--path-repo`.
- `workflow/Snakefile:296` defines the inactive monolithic
  `runSingleAncestryPRS` rule.
- `workflow/scripts/prepare_prs_inputs.sh:31-33` defaults to the vendored tree;
  `:291` and `:296` write `prs_pipeline`/`path_repo`.
- `workflow/scripts/download_prs_resources.sh:132-145` clones/downloads
  software that is already in the image.
- `envs/prs.def:54` bakes the vendored tree into
  `/opt/gdcgenomicsqc/prs_pipeline`.
- `envs/prs.def:76-80` already bakes pinned PRScs, PRScsx, CTSLEB, PROSPER, and
  SDPRX under `/opt/gdcgenomicsqc/software/`.
- `envs/prs.yml:26-28` already installs `prsice`, `plink`, and `plink2`.
- `config/ctsleb_sim_prs.yaml:88-100` overrides image tools with host paths, so
  the image copies are currently not used by that config.

## 5. Implementation steps

### 5.1 Move scripts and make wrappers sibling-relative

1. Move the 12 files listed in section 3.
2. In `prsPipelines.smk`, replace:
   ```python
   PRS_SRC = PRS_PIPELINE_PATH / "src"
   ```
   with:
   ```python
   PRS_SRC = Path(workflow.basedir) / "scripts"
   ```
3. Add a script-relative resolver to shell wrappers:
   ```bash
   SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
   ```
4. Replace hardcoded sibling calls:
   - `workflow/scripts/run_CT.sh:131`
     `Rscript ${path_prs_pipeline}/src/CT_script.R` →
     `Rscript "${SCRIPT_DIR}/CT_script.R"`
   - `workflow/scripts/run_PRScs.sh:150`
     `${prs_pipeline}/src/PRS_sscore_to_R2.R` →
     `${SCRIPT_DIR}/PRS_sscore_to_R2.R`
   - `workflow/scripts/run_PRScsx.sh:192,198`
     `${prs_pipeline}/src/PRS_sscore_to_R2.R` →
     `${SCRIPT_DIR}/PRS_sscore_to_R2.R`
   - `workflow/scripts/score_test.sh:151`
     `${path_repo}/src/evaluate_test.R` →
     `${SCRIPT_DIR}/evaluate_test.R`
5. Remove `--path-repo` parsing and the required `path_repo` check from
   `score_test.sh`; update the `scoreTestPRS` shell block accordingly.
6. Remove the unused `path_repo` positional argument from `run_PRSice2.sh` and
   its rule call.

### 5.2 Remove `prs_pipeline` path plumbing

1. Remove `PRS_PIPELINE_PATH`/`PRS_PIPELINE_PIN` and `checkPRSPipelinePath` from
   `preparePRSInputs.smk`.
2. Remove `--prs-pipeline-dir`, `--prs-pipeline-ref`, `--prs-pipeline-sif`, and
   `--prs-helper-sif` plumbing from `preparePRSInputs.smk`.
3. Remove `path_prs_pipeline` and `prs_pipeline` writes from generated configs in
   `prsPipelines.smk` and `prepare_prs_inputs.sh`.
4. Remove `--prs-pipeline-dir` and the `PRS_PIPELINE_DIR` default from
   `prepare_prs_inputs.sh`.
5. Remove `prs_pipeline_path`, `prs_pipeline_ref`, `single_ancestry_script`, and
   `single_ancestry_flags` from:
   - `config/config.schema.yaml`
   - `config/example_config.yaml`
   - `config/ctsleb_sim_prs.yaml`
   - other tracked examples/smoke configs
6. Delete the inactive monolithic rule `runSingleAncestryPRS`, or replace it with
   a meta-rule that depends on the five native single-method `.done` files.
7. Delete `workflow/scripts/prs_pipeline/run_single_ancestry_PRS_pipeline.sh`.
8. Delete the rest of `workflow/scripts/prs_pipeline/` once the active graph no
   longer references it.

### 5.3 Make image tools authoritative

Change defaults and remove host-path overrides:

| Tool | Target |
|---|---|
| PRSice | `PRSice` from the image’s conda env |
| plink | `plink` from image `PATH` |
| plink2 | `plink2` from image `PATH` |
| PRScs | `/opt/gdcgenomicsqc/software/PRScs` |
| PRScsx | `/opt/gdcgenomicsqc/software/PRScsx` |
| CT-SLEB | `/opt/gdcgenomicsqc/CTSLEB/R` |
| PROSPER | `/opt/gdcgenomicsqc/software/PROSPER` |
| SDPRX | `/opt/gdcgenomicsqc/software/SDPRX` |

Required changes:

1. `runSingleAncestryPRSice`
   - Remove `prsice_r`, `prsice_bin`, and `PRSICE_CMD` assignment.
   - Let `run_PRSice2.sh` use `PRSICE_CMD=${PRSICE_CMD:-PRSice}`.
2. `scoreTestPRS`
   - Remove the `prsice_bin` staging symlink.
   - Ensure `PRSice` is on `PATH` from the image.
3. `single_prscs` / `multi_prscsx`
   - Default `path_code` to `/opt/gdcgenomicsqc/software/PRScsx`.
   - Keep `ld_ref_dir` configurable as read-only data.
4. `multi_ctsleb`
   - Default `software_dir` to `/opt/gdcgenomicsqc/CTSLEB/R`.
5. `multi_prosper`
   - Default `software_dir` to `/opt/gdcgenomicsqc/software/PROSPER`.
   - Default `plink2` to `plink2` on the image `PATH`.
   - Add `plink2` to the `singleprshelper` image, or switch the PROSPER rules to
     the `prs` image. Do not fall back to
     `/projects/standard/gdc/public/plink2`.
6. Remove host-path overrides from `config/ctsleb_sim_prs.yaml` so image paths are
   actually exercised.

### 5.4 Update `download_prs_resources.sh` and resource checks

1. Remove software cloning/downloading:
   - PRScs, PRScsx, CTSLEB, PROSPER, SDPRX clones.
   - PRSice archive download.
   - `--download-software`.
2. Remove `prs_pipeline` verification and provenance recording.
3. Keep only reference-data provisioning if LD panels remain external.
4. Replace `preparePRSMethodResources` with a lightweight
   `checkPRSImageTools` rule, or keep the rule name but make its script verify
   image tools:
   - `command -v PRSice`
   - `command -v plink`
   - `command -v plink2`
   - PRScsx code path
   - CTSLEB path
   - PROSPER path
5. Replace `resources.ready` with a marker such as
   `OUT_DIR/.prs_image_tools.checked`, recording image paths and versions.
6. Update every rule that currently depends on
   `rules.preparePRSMethodResources.output.ready`.

### 5.5 Update the PRS image

In `envs/prs.def`:

1. Remove:
   ```dockerfile
   ../workflow/scripts/prs_pipeline /opt/gdcgenomicsqc/prs_pipeline
   ```
2. Keep:
   ```dockerfile
   ../workflow/scripts/CTSLEB /opt/gdcgenomicsqc/CTSLEB
   ```
3. Keep the pinned clones under `/opt/gdcgenomicsqc/software/`.
4. Optionally bake the 12 flat scripts under
   `/opt/gdcgenomicsqc/scripts/` for image-only execution. The canonical copies
   remain in `workflow/scripts/`.
5. Update the image comments and `envs/README.md`.
6. Rebuild and push only after the local smoke/end-to-end path passes.

### 5.6 Update configs, schema, and smoke harness

1. Remove host tool paths from:
   - `config/example_config.yaml`
   - `config/ctsleb_sim_prs.yaml`
   - tracked MSI/sandbox examples
2. Update schema descriptions for `path_code`, `software_dir`, `plink2`, and
   `resource_dir`.
3. Update `workflow/scripts/smoke_single_ancestry.sh:109` so it checks the flat
   scripts instead of `prs_pipeline/run_single_ancestry_PRS_pipeline.sh`.
4. Update `AGENTS.md`, `docs/PRS_full_integration_plan.md`, and `envs/README.md`
   to describe flat scripts and image-provided tools.

## 6. Verification

Run all commands from `workflow/` unless noted otherwise.

### 6.1 Static checks

```bash
rg -n "prs_pipeline|PRS_PIPELINE|path_repo|PRSice\\.R|PRSice_linux" \
  workflow/rules workflow/scripts workflow/Snakefile config envs docs
```

Expected: no active references remain.

### 6.2 Image checks

Inside `envs/prs.sif`:

```bash
command -v PRSice
command -v plink
command -v plink2
test -f /opt/gdcgenomicsqc/software/PRScsx/PRScsx.py
test -d /opt/gdcgenomicsqc/CTSLEB/R
test -d /opt/gdcgenomicsqc/software/PROSPER/scripts
```

### 6.3 Snakemake dry runs

```bash
/projects/standard/gdc/public/envs/snakemake/bin/snakemake -n run_allPRSPipelines \
  --rerun-incomplete \
  --profile ../profiles/hpc \
  --configfile ../config/ctsleb_sim_prs.yaml

/projects/standard/gdc/public/envs/snakemake/bin/snakemake -n run_scoreTestPRS \
  --rerun-incomplete \
  --profile ../profiles/hpc \
  --configfile ../config/ctsleb_sim_prs.yaml
```

### 6.4 Smoke workflow

```bash
sbatch workflow/scripts/smoke_single_ancestry.sh
```

### 6.5 Full end-to-end run

```bash
/projects/standard/gdc/public/envs/snakemake/bin/snakemake run_allPRSPipelines \
  --rerun-incomplete \
  --keep-going \
  --profile ../profiles/hpc \
  --configfile ../config/ctsleb_sim_prs.yaml
```

Launch this under a durable allocation (`sbatch`/`tmux`); do not run it inside a
1-day OnDemand interactive job. Complete the fixes in section 6.6 first.

Success criteria:

- All 10 `method_runs/*.done` markers exist.
- `test_evaluation/` contains all expected results/scores.
- Output file names and R² values match the pre-move run.
- No `prs_pipeline` path appears in logs or generated configs.

### 6.6 Full-scale `ctsleb_sim_prs` remediation (2026-11-02 run)

The 2026-11-02 full-scale run did not complete. It was not a de-vendoring
regression; it exposed three pre-existing full-scale issues plus two operational
issues. Resolve these before repeating section 6.5.

#### 6.6.1 Observed state

- The rerun command used `run_allPRS_Pipelines` (with an underscore) instead of
  the real target `run_allPRSPipelines` (`workflow/Snakefile:349`). No new
  Snakemake log was created for 2026-11-03, so the intended rerun never started.
- The 2026-11-02 orchestrator was running inside a 1-day OnDemand interactive
  allocation (`5546127`, TIMEOUT). The main Snakemake log ends at 13:37 while
  child SLURM jobs kept running unsupervised.
- Completed jobs: PRSice (`5549508`, 3h03), single-ancestry PRS-CS (`5549504`,
  9h20), PROSPER prep (`5549502`, 45m).
- Failed/timed-out jobs: CT (`5549510`, TIMEOUT at 4h), LD matrix (`5549505`,
  TIMEOUT at 10h), PRS-CSx (`5549501`, TIMEOUT at 12h), CTSLEB (`5549507`,
  FAILED after 7h02).
- `.done` markers present: `single_prsice.done`, `single_prscs.done`,
  `multi_prosper.done`, plus stale adapter stubs `multi_ldpred2.done` and
  `multi_sdprs.done` from 2026-10-26.
- Missing: `single_ct.done`, `single_ldpred2.done`, `single_lassosum2.done`,
  `multi_ctsleb.done`, `multi_prscsx.done`, and `test_evaluation/`.

#### 6.6.2 Root causes

1. **Wrong target:** `run_allPRS_Pipelines` does not exist.
2. **Orchestrator died:** a 1-day interactive allocation expired; no scheduler
   remained to continue failed or pending jobs.
3. **CTSLEB `.sscore` column mismatch:** `helper_CombinePRS.R:25` selects
   `prs_temp[, 5:ncol(prs_temp)]`, but current plink2
   (`v2.0.0-a.6.12LM`) emits `#FID IID PHENO1 ALLELE_CT
   NAMED_ALLELE_DOSAGE_SUM SCORE1_SUM ... SCORE12_SUM` (17 columns). The score
   columns start at 6, not 5. The one-column shift pads names with `NA`,
   `bind_cols` repairs them to `...13`, `...26`, etc., and `GetSNPSet` fails on
   a name with no underscores (`str_temp[[1]][[8]]`: subscript out of bounds).
4. **LD matrix near-dense:** `generate_ld_matrix.R:30` uses
   `map$genetic.dist`, which is all zeros in the CTSLEB `.bim` files, with
   `size = 3/1000`. `snp_cor` therefore treats every SNP on a chromosome as
   colocated, producing a chromosome-wide matrix instead of a 3 cM sparse
   matrix. After 10h, no `chr*_corr.rds` had been written. The preferred
   workaround is to supply a precomputed external matrix through
   `prsMethods.<method>.ld_matrix_dir` (section 6.6.3 F); the window fallback
   remains the in-pipeline fix.
5. **Runtimes too short for full scale:** CT was still in `--indep-pairwise` at
   334k/1.9M markers after 4h (needs ~23h+); PRS-CSx reached chr12 after 12h
   (needs ~22h+); CTSLEB had a 12h limit.

#### 6.6.3 Code and resource fixes

Applied to the active copies on 2026-11-03. The de-vendoring move in section 5.1
must carry the same fixes into `workflow/scripts/generate_ld_matrix.R`,
`workflow/scripts/run_LDpred2.R`, and `workflow/scripts/run_lassosum2.R`.

**A. CTSLEB score-column selection** — `workflow/scripts/CTSLEB/R/helper_CombinePRS.R:25`

```r
score_cols <- grep("^SCORE[0-9]+_SUM$", colnames(prs_temp))
if (length(score_cols) == 0L) {
  score_cols <- 5:ncol(prs_temp)
  warning("No SCORE*_SUM columns found; falling back to columns 5:ncol")
}
if (length(score_cols) != length(names)) {
  stop(sprintf("Expected %d score columns, found %d",
               length(names), length(score_cols)))
}
prs_list[[temp]] <- prs_temp[, score_cols, drop = FALSE]
```

Also add a guard in `workflow/scripts/run_ctsleb_wrapper.R:123` after
`best_snps` is selected:

```r
if (!grepl("^clump_r2_.+_ws_.+_p_other_.+_p_tar_", best_snps)) {
  stop("CTSLEB wrapper: unexpected best_snps name: ", best_snps)
}
```

**B. LD matrix window** — `workflow/scripts/prs_pipeline/src/generate_ld_matrix.R:30,46`

```r
gen <- obj.bigSNP$map$genetic.dist
use_genetic <- !is.null(gen) && !all(is.na(gen)) && any(gen != 0, na.rm = TRUE)
POS2 <- if (use_genetic) gen else obj.bigSNP$map$physical.pos
window <- if (use_genetic) 3 / 1000 else 3000  # 3 cM or 3 Mb
...
corr0 <- snp_cor(G, ind.col = ind.chr, size = window,
                 infos.pos = POS2[ind.chr], ncores = args$ncores)
```

The same `use_genetic`/`window` fallback is applied to the internal `snp_cor`
branches in `run_LDpred2.R:212-224` and `run_lassosum2.R:237-244` (the pipeline
normally uses the precomputed matrix, but these branches carried the same latent
bug).

**C. Resource raises** — `workflow/rules/prsPipelines.smk`,
`workflow/rules/ctsleb_multi_prs_rule.smk:28`

| Rule | Change |
|---|---|
| `generateLDMatrix` | `mem_mb=64000`, `runtime=1440` |
| `runSingleAncestryCT` | `mem_mb=32000`, `runtime=2880` |
| `runSingleAncestryLDpred2` | `runtime=2160` |
| `runSingleAncestryLassosum2` | `runtime=2160` |
| `runMultiAncestryPRSCSx` | `runtime=1800` |
| `runMultiAncestryCTSLEB` | `runtime=1440` |

**D. External test sample with `test_fraction: 0`** —
`workflow/scripts/prepare_prs_inputs.sh`

Allow `--test-fraction 0` when `--test-prefix` is supplied. Only carve a test set
and remove it from ANC1 when the fraction is greater than zero; otherwise use
the supplied test prefix directly. This is the path the toy smoke harness uses.

**E. Durable orchestrator** — write a small `sbatch` wrapper (or use
`sbatch --wrap` with `--chdir=.../workflow`) on `msilong` with a long wall time,
and keep the orchestrator out of a 1-day interactive allocation.

**F. External precomputed LD matrix (preferred for full-scale)** — upstream
`prs_pipeline/docs/sandbox_singularity_runner.md` documents `ld_matrix_dir` as a
precomputed directory produced by `src/generate_ld_matrix.R`; it takes priority
over `ld_cache_dir`. The GDC wrapper exposes the same directory per method, but
the full/reduced configs did not set it, so `generateLDMatrix` ran from the
simulated genotypes and timed out.

Set the per-method override (not the shared `prsPipeline.ld_matrix_dir`, which is
the `generateLDMatrix` output directory). Add these entries under the existing
`prsMethods` block; do not duplicate the key:

```yaml
prsMethods:
    single_ldpred2:
        ld_matrix_dir: "/projects/standard/gdc/public/prs_methods/ref/ldpred2"
    single_lassosum2:
        ld_matrix_dir: "/projects/standard/gdc/public/prs_methods/ref/ldpred2"
```

The directory must contain:

- `map.rds` with `chr`, `rsid`, `pos`, `a1`, `a0` (optionally `af_UKBB`, `ld`)
  in the same row order as the concatenated chromosome matrices.
- `chr1_corr.rds` ... `chr22_corr.rds`. LDpred2 reads only this filename;
  lassosum2 also accepts `chr<N>.rds` / `LD_chr<N>.rds`.
- `g_idx.rds` only when the directory is also used as the shared
  `prsPipeline.ld_matrix_dir`. LDpred2/lassosum2 do not read `g_idx.rds`; using
  per-method overrides avoids pulling `generateLDMatrix` when it is absent.

The matrix may come from an external reference panel (here the LDpred2/UKBB HM3
matrix under `ref/ldpred2`), not the study genotypes. LDpred2 matches summary
stats to `map.rds` by RSID and then restricts to variants present in the target
genotypes. For the CTSLEB simulation, about 1.02M of 1.9M variants overlap this
HM3 reference, which is expected. Verify each chromosome matrix's row count
against `map.rds` before launching.

To generate a matrix outside the pipeline:

```bash
Rscript workflow/scripts/prs_pipeline/src/generate_ld_matrix.R \
  --anc_bed /path/to/study.bed \
  --out /path/to/external_ld_matrix \
  --ncores 16
```

The GDC rules pass only `--ld-matrix-dir`; `ld_cache_dir` remains an upstream
option and is not wired into the Snakemake rules. Keep `afreq_file` configured
for LDpred2/lassosum2; it is independent of the LD matrix. `ref/ldpred2` lives
under `/projects`, which the `hpc` profile already binds; add other roots with
`--bind`. Because the external directory is passed as a parameter rather than a
declared input, changing its contents does not automatically invalidate the
method outputs; remove the affected `.done` markers or use `--forcerun`.

#### 6.6.4 Phase 1 — reduced end-to-end validation

Run the toy smoke in section 6.6.7 first; it is the fastest regression check for
the code fixes. Then create a local, gitignored
`config/ctsleb_sim_prs_small.yaml` from `config/ctsleb_sim_prs.yaml` with only
these changes:

```yaml
OUT_DIR: "/scratch.global/baron063/testing/ctsleb_sim_prs_small"
local-storage-prefix: "/scratch.global/baron063/testing/ctsleb_sim_prs_small/tmp"

phenotypeSimulation:
    skip_thinning: false
    thin_count_snps: 200000
    thin_count_inds: 10000

prsPipeline:
    n_total_gwas: 4000
    test_fraction: 0.2
    test_sample: "/scratch.global/baron063/testing/ctsleb_sim_prs_small/prs_inputs/AFR_EUR/test_sample/AFR_simulation_test_sample"
    generated_input_dir: "/scratch.global/baron063/testing/ctsleb_sim_prs_small/prs_inputs/AFR_EUR"

prsMethods:
    single_ldpred2:
        ld_matrix_dir: "/projects/standard/gdc/public/prs_methods/ref/ldpred2"
    single_lassosum2:
        ld_matrix_dir: "/projects/standard/gdc/public/prs_methods/ref/ldpred2"
```

`n_total_gwas: 4000` reflects the actual split: 10k simulated individuals → 2k
test → 8k remaining → 4k GWAS / 4k study. Keep the original bfile
`input_prefixes`; do not point them at the existing full-scale `.pgen` files,
because `runPhenotypeSimulation.R` overwrites its input prefix when thinning.

Run from `workflow/` under `sbatch`/`tmux`:

```bash
/projects/standard/gdc/public/envs/snakemake/bin/snakemake -n run_allPRSPipelines \
  --profile ../profiles/hpc \
  --configfile ../config/ctsleb_sim_prs_small.yaml

/projects/standard/gdc/public/envs/snakemake/bin/snakemake run_allPRSPipelines \
  --rerun-incomplete \
  --keep-going \
  --profile ../profiles/hpc \
  --configfile ../config/ctsleb_sim_prs_small.yaml
```

Expected wall time: ~2–4h, with tens of GB of temporary space for the
bfile→pfile conversion. This run is the regression check for the CTSLEB fix, the
external-matrix wiring, the reduced-data path, and the full method graph; the
internal LD-matrix window fallback is covered by the toy smoke.

#### 6.6.5 Phase 2 — full-scale rerun

Reuse `config/ctsleb_sim_prs.yaml` and the existing OUT_DIR. Add the per-method
`ld_matrix_dir` block from section 6.6.3 F (it is already present in the local
config) so the run uses the existing `ref/ldpred2` matrix instead of regenerating
it. After Phase 1 is green:

```bash
/projects/standard/gdc/public/envs/snakemake/bin/snakemake -n run_allPRSPipelines \
  --profile ../profiles/hpc \
  --configfile ../config/ctsleb_sim_prs.yaml
```

Then launch the durable orchestrator with:

```bash
/projects/standard/gdc/public/envs/snakemake/bin/snakemake run_allPRSPipelines \
  --rerun-incomplete \
  --keep-going \
  --profile ../profiles/hpc \
  --configfile ../config/ctsleb_sim_prs.yaml
```

With the external matrix, the dry run should **not** contain
`generateLDMatrix`; it should start at `preparePRSInputs`/the method rules. If
`generateLDMatrix` still appears, one of the two methods is still using the
shared/default directory. Verified 2026-11-03: the full config now builds 7
jobs (`single_ct`, `single_ldpred2`, `single_lassosum2`, `multi_ctsleb`,
`multi_prscsx`, `scoreTestPRS`, `run_allPRSPipelines`), and the reduced config
builds 25 jobs including PROSPER and the adapter stubs; neither dry run contains
`generateLDMatrix`. `--rerun-incomplete` should rerun only `single_ct`,
`single_ldpred2`, `single_lassosum2`, `multi_ctsleb`, `multi_prscsx`, and
`scoreTestPRS`. `single_prsice`, `single_prscs`, and `multi_prosper` are already
complete and should be skipped. Expected wall time: ~1.5–2.5 days, dominated by
CT and PRS-CSx.

#### 6.6.6 Verification

- Exactly 10 `method_runs/*.done` markers exist.
- `test_evaluation/` contains the expected scores/results.
- CTSLEB `*_best_snps.tsv` contains a valid
  `clump_r2_..._ws_..._p_other_..._p_tar_...` name, and no `...N` columns appear
  in `prs_mat`.
- With an external matrix, `generateLDMatrix` is absent from the DAG; the
  `runSingleAncestryLDpred2.log` and `runSingleAncestryLassosum2.log` files start
  with `Loading pre-computed LD matrix from: ...` and retain more than 10
  variants after matching.
- If the internal generator is used instead, it writes `map.rds`, `g_idx.rds`,
  and all `chr*_corr.rds`; `length(corr@x) / length(corr)` is below the
  chromosome-wide density, not equal to it.
- CT finishes within 48h; PRS-CSx reaches chr22; CTSLEB no longer raises
  `GetSNPSet`/`str_temp[[1]][[8]]`.
- No `prs_pipeline` path appears in logs or generated configs after the
  de-vendoring changes.

#### 6.6.7 Toy smoke validation (2026-11-03)

The toy smoke passed after the section 6.6.3 fixes. It uses 300 real chr22
HapMap3 SNPs and 100 samples per ancestry, with 50 study samples after the
50/50 GWAS/study split, plus an external held-out test sample. Wall time was
7m05s on 8 cores and 32 GB.

```bash
sbatch --chdir="$PWD" --cpus-per-task=8 --mem=32G --time=02:00:00 \
  --wrap='bash workflow/scripts/smoke_single_ancestry.sh'
```

Use `--wrap` (or `srun`) rather than `sbatch workflow/scripts/smoke_single_ancestry.sh`:
Slurm copies the latter to `/var/spool`, so `BASH_SOURCE` no longer points at the
repo and the harness cannot find `workflow/Snakefile`.

Results:

- All five single-ancestry methods, `scoreTestPRS`, `multi_prscsx`, and
  `multi_ctsleb` passed.
- CTSLEB produced `clump_r2_0.01_ws_5000_p_other_5e-08_p_tar_5e-08` and did not
  raise `GetSNPSet`/`str_temp[[1]][[8]]`.
- `generateLDMatrix` wrote `chr22_corr.rds` (300 SNPs, density 0.499) without
  stalling. The 3 Mb physical window is relatively dense on a 300-SNP toy set
  because the SNPs are concentrated in a small region, but it is no longer
  chromosome-wide.
- The external `test_sample` + `test_fraction: 0` path worked after the
  `prepare_prs_inputs.sh` fix.

Not covered by the toy smoke: PROSPER, the `multi_ldpred2`/`multi_sdprs` adapter
stubs, and full-scale runtime behavior. Those remain for sections 6.6.4/6.6.5.

## 7. Estimate

| Work item | Estimate |
|---|---|
| Move 12 scripts and rewire sibling paths | 0.5–1 day |
| Remove path plumbing and monolithic rule | 0.5–1 day |
| Image tool defaults, PRSice switch, PROSPER image fix | 1 day |
| Remove software downloads, update resource checks | 0.5 day |
| Image rebuild/push | 0.5 day |
| Dry runs, smoke, and full end-to-end validation | 1–2 days |
| Full-scale `ctsleb_sim_prs` remediation (section 6.6) | 0.5 day code/resources |

**Total: approximately 3–4 focused days**, assuming LD references stay externally
mounted. Add 1–3 days if LD references must also be baked into the image. Section
6.6 adds ~0.5 day of code/resource work plus ~2–3 days of wall-clock time for the
reduced and full-scale validation runs.

## 8. Risks and mitigations

| Risk | Mitigation |
|---|---|
| Image PRSice version differs from bundled PRSice.R behavior | Run PRSice on the smoke dataset; compare `.prsice`, `.best`, and `.snps` outputs |
| `singleprshelper` lacks `plink2` | Add plink2 to that image or switch PROSPER to the `prs` image |
| Stale generated configs still contain `prs_pipeline` | Grep `OUT_DIR` and generated `*.conf` before full run; delete stale configs or rerun `preparePRSInputs` |
| Image path drift | Centralize image tool paths in `prsPipelines.smk` constants and verify in one marker rule |
| Script relocation breaks relative calls | Use `SCRIPT_DIR` in every shell wrapper; test each method individually |
| Existing configs override image paths | Remove host overrides from tracked configs; keep they only as optional legacy overrides |
| plink2 `.sscore` baseline column layout changes | Select score columns by name (`^SCORE[0-9]+_SUM$`), never by fixed indices; validate the count in the CTSLEB helper |
| PLINK `.bim` genetic distances are all zero | Fall back to physical positions and a bp-based window (3 Mb) in every `snp_cor` call |
| Full-scale methods exceed default SLURM runtimes | Apply the resource table in section 6.6.3 and run under a durable orchestrator, not a 1-day interactive allocation |

## 9. Rollback

1. Keep the vendored tree until the flat scripts pass smoke and end-to-end tests.
2. Make the move in two commits:
   - Add flat copies + rewire paths; validate.
   - Delete vendored tree and path plumbing; validate again.
3. Revert the second commit alone if cleanup causes issues.

## 10. Open questions

1. Should LD reference panels remain externally mounted data, or be baked into the
   singularity images?
2. Should `run_singleAncestryPRS` remain as a compatibility meta-target, or be
   removed entirely?
3. Is the `prs` image sufficient for PROSPER, or must `singleprshelper` be updated
   with `plink2`?
