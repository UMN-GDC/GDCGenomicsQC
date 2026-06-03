# Lab: Visualizing and Working with Single-Ancestry PRS Outputs
#
# Interact with output files from the Single-Ancestry PRS Pipeline.
# Work with PRS scores, performance metrics (R2, AUC), and compare methods.
#
# Prerequisites:
# - Completed the PRS tutorial (outputs in OUT_DIR/prs_inputs/{ANC1}_{ANC2}/)
# - R (>= 4.0) with tidyverse installed
# - Set OUT_DIR to your pipeline output directory

conda_env <- '/scratch.global/R25_files/r25_r_env/lib/R/library'
.libPaths(unique(c(.Library.site, .Library)))
.libPaths(c(conda_env, .libPaths()))

library(tidyverse)

OUT_DIR <- "/path/to/your/pipeline/output"
anc1 <- "AFR"
anc2 <- "EUR"
prs_out_dir <- file.path(OUT_DIR, "prs_inputs", paste0(anc1, "_", anc2))

# Load PRS scores for PRSice2
prsice_score_path <- file.path(prs_out_dir, "single_prsice", "best_pRS.prs")
if (file.exists(prsice_score_path)) {
  prs_scores <- read_tsv(prsice_score_path, col_names = c("FID", "IID", "PRS"))
  glimpse(prs_scores)
}

# Load PRSice2 summary
prsice_summary_path <- file.path(prs_out_dir, "single_prsice", "summary.csv")
if (file.exists(prsice_summary_path)) {
  prsice_summary <- read_csv(prsice_summary_path)
  glimpse(prsice_summary)
}

# Load LDPred2 scores
ldpred2_path <- file.path(prs_out_dir, "single_ldpred2", "prs_scores.tsv")
if (file.exists(ldpred2_path)) {
  ldpred2_scores <- read_tsv(ldpred2_path)
  glimpse(ldpred2_scores)
}

# PRS score distribution (PRSice2)
if (exists("prs_scores")) {
  ggplot(prs_scores, aes(x = PRS)) +
    geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
    labs(
      title = "PRS Score Distribution (PRSice2)",
      x = "Polygenic Risk Score",
      y = "Number of Samples"
    ) +
    theme_minimal()
}
