#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(caret)
  library(SuperLearner)
  library(ranger)
  library(glmnet)
  library(pROC)
})

parse_cli_args <- function(argv, required) {
  if (length(argv) %% 2 != 0L || any(!startsWith(argv[c(TRUE, FALSE)], "--"))) {
    stop("Arguments must be supplied as --name value pairs.")
  }

  keys <- gsub("-", "_", sub("^--", "", argv[c(TRUE, FALSE)]))
  values <- argv[c(FALSE, TRUE)]
  if (anyDuplicated(keys)) {
    stop("Duplicate command-line argument(s): ", paste(unique(keys[duplicated(keys)]), collapse = ", "))
  }
  missing <- setdiff(required, keys)
  if (length(missing) > 0L) {
    stop("Missing required argument(s): --", paste(gsub("_", "-", missing), collapse = ", --"))
  }
  as.list(stats::setNames(values, keys))
}

ctsleb_src <- "/scratch.global/saonli/GDCGenomicsQC/prs_resources/software/CTSLEB_github/R"
r_files <- list.files(ctsleb_src, pattern = "\\.[Rr]$", full.names = TRUE)

if (length(r_files) == 0) {
  stop("No CTSLEB R source files found in: ", ctsleb_src)
}

invisible(lapply(sort(r_files), source))

args <- parse_cli_args(
  commandArgs(trailingOnly = TRUE),
  required = c(
    "target_sumstats", "aux_sumstats", "target_ref_plink", "aux_ref_plink",
    "target_test_plink", "plink19", "plink2", "results_dir", "out_prefix",
    "tuning_pheno", "validation_pheno"
  )
)

dir.create(args$results_dir, recursive = TRUE, showWarnings = FALSE)

sum_ref <- fread(args$aux_sumstats)
sum_target <- fread(args$target_sumstats)

params_farm <- SetParamsFarm(
  plink19_exec = args$plink19,
  plink2_exec = args$plink2
)

# Step 1: two-dimensional clumping and thresholding.
prs_mat <- dimCT(
  results_dir = args$results_dir,
  sum_target = sum_target,
  sum_ref = sum_ref,
  ref_plink = args$aux_ref_plink,
  target_plink = args$target_ref_plink,
  test_target_plink = args$target_test_plink,
  out_prefix = args$out_prefix,
  params_farm = params_farm
)

y_tune <- fread(args$tuning_pheno)
y_valid <- fread(args$validation_pheno)

n_tune <- nrow(y_tune)
if (nrow(prs_mat) < (n_tune + nrow(y_valid))) {
  stop("CTSLEB wrapper: PRS matrix rows are fewer than tuning + validation samples.")
}

prs_tune <- prs_mat[seq_len(n_tune), ]
prs_valid <- prs_mat[seq.int(n_tune + 1, n_tune + nrow(y_valid)), ]

# The tutorial computes tuning R^2 across all PRS columns.
prs_cols <- seq.int(3, ncol(prs_tune))
prs_r2 <- vapply(
  prs_cols,
  function(i) summary(lm(y_tune[[1]] ~ prs_tune[[i]]))$r.squared,
  numeric(1)
)
max_ind <- which.max(prs_r2)
best_snps <- colnames(prs_tune)[max_ind + 2]

# Step 2: empirical Bayes effect sizes.
prs_mat_eb <- CalculateEBEffectSize(
  bfile = args$target_test_plink,
  snp_ind = best_snps,
  plink_list = plink_list,
  out_prefix = args$out_prefix,
  results_dir = args$results_dir,
  params_farm = params_farm
)

# Step 3: clean PRS matrices and store intermediates for a later SL step.
prs_tune_eb <- prs_mat_eb[seq_len(n_tune), ]
prs_valid_eb <- prs_mat_eb[seq.int(n_tune + 1, n_tune + nrow(y_valid)), ]

# Lock the empirical-Bayes model using tuning outcomes only, then report the
# corresponding validation performance once.  The first two columns are the
# participant identifiers, matching the dimCT output convention.
eb_cols <- seq.int(3, ncol(prs_tune_eb))
eb_cols <- eb_cols[vapply(eb_cols, function(i) {
  values <- suppressWarnings(as.numeric(prs_tune_eb[[i]]))
  all(is.finite(values)) && stats::sd(values) > 0
}, logical(1))]
if (length(eb_cols) == 0L) {
  stop("CTSLEB wrapper: no finite, non-constant empirical-Bayes PRS columns were produced.")
}

y_tune_values <- suppressWarnings(as.numeric(y_tune[[1]]))
y_valid_values <- suppressWarnings(as.numeric(y_valid[[1]]))
if (any(!is.finite(y_tune_values)) || any(!is.finite(y_valid_values))) {
  stop("CTSLEB wrapper: tuning or validation phenotype contains non-finite values.")
}

eb_tuning_r2 <- vapply(eb_cols, function(i) {
  summary(lm(y_tune_values ~ as.numeric(prs_tune_eb[[i]])))$r.squared
}, numeric(1))
best_eb_col <- eb_cols[which.max(eb_tuning_r2)]
best_eb_name <- colnames(prs_tune_eb)[best_eb_col]
validation_model <- lm(y_valid_values ~ as.numeric(prs_valid_eb[[best_eb_col]]))
validation_summary <- summary(validation_model)
validation_coef <- validation_summary$coefficients[2, ]
fair_metrics <- data.table(
  method = "CT-SLEB",
  selected_prs_column = best_eb_name,
  tuning_r2 = max(eb_tuning_r2),
  test_r2 = validation_summary$r.squared,
  test_beta = validation_coef[["Estimate"]],
  test_se = validation_coef[["Std. Error"]],
  test_p = validation_coef[["Pr(>|t|)"]],
  test_n = nrow(prs_valid_eb)
)
fair_prs <- data.table(
  FID = prs_valid_eb[[1]],
  IID = prs_valid_eb[[2]],
  PRS = as.numeric(prs_valid_eb[[best_eb_col]])
)

cleaned <- PRS_Clean(
  Tune_PRS = prs_tune_eb,
  Tune_Y = y_tune,
  Validation_PRS = prs_valid_eb
)

fwrite(as.data.table(prs_mat), file.path(args$results_dir, paste0(args$out_prefix, "_dimct_prs.tsv")), sep = "\t")
fwrite(as.data.table(prs_mat_eb), file.path(args$results_dir, paste0(args$out_prefix, "_eb_prs.tsv")), sep = "\t")
fwrite(as.data.table(cleaned$Cleaned_Tune_PRS), file.path(args$results_dir, paste0(args$out_prefix, "_tune_clean.tsv")), sep = "\t")
fwrite(as.data.table(cleaned$Cleaned_Validation_PRS), file.path(args$results_dir, paste0(args$out_prefix, "_validation_clean.tsv")), sep = "\t")
fwrite(data.table(best_snps = best_snps), file.path(args$results_dir, paste0(args$out_prefix, "_best_snps.tsv")), sep = "\t")
fwrite(fair_metrics, file.path(args$results_dir, paste0(args$out_prefix, "_fair_metrics.tsv")), sep = "\t")
fwrite(fair_prs, file.path(args$results_dir, paste0(args$out_prefix, "_test_prs.tsv")), sep = "\t")

cat("CTSLEB draft wrapper completed.\n")
cat("Best SNP set:", best_snps, "\n")
cat("Selected EB PRS:", best_eb_name, "\n")
cat("Held-out validation R2:", fair_metrics$test_r2, "\n")
