# CTSLEB vendored sources

These 30 `R/*.R` source files are the CT-SLEB multi-ancestry PRS method
(2D clumping/thresholding + empirical Bayes + super-learning), vendored so
the workflow does not depend on a hardcoded scratch path.

- Upstream: andrewhaoyu/CTSLEB (R package, LICENSE: MIT, Copyright (c) 2022
  Haoyu Zhang; see `LICENSE`).
- Source of the copy: the CTSLEB package checkout at
  `/scratch.global/saonli/GDCGenomicsQC/prs_resources/software/CTSLEB_github`
  (as mirrored on SARS-HPC). Files are byte-identical to that checkout.
- Date vendored: 2026-09-28.

The `prs` image also carries a copy at `/opt/gdcgenomicsqc/CTSLEB/R`; the
workflow defaults to this repo copy because it includes the plink2 `.sscore`
`SCORE*_SUM` fix (docs/PRS_de_vendoring_plan.md §6.6.3 A). Point
`prsMethods.multi_ctsleb.software_dir` at the image copy only after rebuilding
the image with that fix.

Main entry points used by `run_ctsleb_wrapper.R`: `SetParamsFarm()`,
`dimCT()`, `CalculateEBEffectSize()`, `PRS_Clean()`, `ExtractFinalBetas()`.
`dimCT()` assigns the globals `sum_com`, `write_list`, `snp_list`,
`plink_list`, and `prs_mat`; `CalculateEBEffectSize()` assigns `best_snps_set`,
`unique_infor_post`, `plink_list_eb`, `scores_eb`, `score_eb_file`,
`p_values_eb`.