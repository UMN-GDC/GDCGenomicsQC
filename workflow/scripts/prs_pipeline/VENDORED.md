# Vendored `prs_pipeline` (absorbed — no longer maintained upstream)

These files are an **absorbed copy** of the `UMN-GDC/prs_pipeline` method engine.
Upstream is treated as frozen: it will not be used or updated from here; any
changes happen in GDCGenomicsQC as native Snakemake rules.

Source of this copy:
- Repo: `git@github.com:UMN-GDC/prs_pipeline.git`
- Branch/commit: `sandbox_multi_pheno` @ `f71cf4f1031ff6231bc9f8abc9f96ebdc6f17bdf`
  (2026-08-20) — same commit recorded as `prs_pipeline_ref` in the workflow config.
- Vendored into GDCGenomicsQC on 2026-09-24 from `git archive` of that commit.

What is vendored (matches the commit exactly):
- `src/` — per-method scripts (`run_CT.sh`, `run_LDpred2.R`, `run_lassosum2.R`,
  `run_PRSice2.sh`, `prepare_sumstats.R`, `generate_ld_matrix.R`, `score_test.sh`,
  PRS-CSx/VIPRS/PROSPER runners, data-prep splitter, `PRSice.R` + `PRSice_linux`)
- `templates/` — shell-readable config templates
- `envs/` — Singularity `.def` + conda `.yaml` for `singlePRS`, `multiPRS`, `singlePRSHelper`
- `run_single_ancestry_PRS_pipeline.sh` — monolithic orchestrator, kept only as a
  transitional fallback for `rule runSingleAncestryPRS` until Phase 2 replaces it
  with native per-method rules.
- `LICENSE` (MIT, UMN-GDC)

Intentionally NOT vendored: `temp_ld/` (runtime LD cache), `prsv2_latest.sif` /
`singleprshelper_latest.sif` (pull from GHCR, see `prsMethods.containers.*`), `docs/`,
README, and the deprecated orchestrators (`run_prepare_prs.sh`, `run_split_data.sh`,
`sandbox_singularity_runner.sh`, …) — these are reimplemented as rules.

Do not hand-edit vendored scripts; modify them via the normal GDCGenomicsQC
workflow/normal branches so the engine stays auditably in lockstep with the repo.