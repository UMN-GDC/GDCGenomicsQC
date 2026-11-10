#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(dplyr)
  library(caret)
  library(SuperLearner)
  library(ranger)
  library(glmnet)
  library(pROC)
})

# Locate the vendored CTSLEB R source directory unless overridden.
resolve_ctsleb_src <- function(path) {
  if (!is.null(path) && path != "" && dir.exists(path)) {
    return(path)
  }
  env_dir <- Sys.getenv("CTSLEB_R_DIR", "")
  if (env_dir != "" && dir.exists(env_dir)) {
    return(env_dir)
  }
  args0 <- commandArgs(trailingOnly = FALSE)
  script <- sub("^--file=", "", args0[grepl("^--file=", args0)][1])
  if (length(script) && script != "" && dir.exists(script)) {
    script <- file.path(script, "run_ctsleb_wrapper.R")
  }
  vendored <- file.path(dirname(normalizePath(script)), "CTSLEB", "R")
  if (dir.exists(vendored)) {
    return(vendored)
  }
  stop("Cannot find CTSLEB R sources (pass --ctsleb-src, set CTSLEB_R_DIR, or use the vendored copy).")
}

parser <- ArgumentParser(description = "CT-SLEB multi-ancestry PRS wrapper for GDCGenomicsQC")
parser$add_argument("--target-sumstats", required = TRUE)
parser$add_argument("--aux-sumstats", required = TRUE)
parser$add_argument("--target-ref-plink", required = TRUE)
parser$add_argument("--aux-ref-plink", required = TRUE)
parser$add_argument("--target-test-plink", required = TRUE)
parser$add_argument("--plink19", required = TRUE)
parser$add_argument("--plink2", required = TRUE)
parser$add_argument("--results-dir", required = TRUE)
parser$add_argument("--out-prefix", required = TRUE)
parser$add_argument("--tuning-pheno", required = TRUE)
parser$add_argument("--validation-pheno", required = TRUE)
parser$add_argument("--ctsleb-src", default = "", help = "directory of the CTSLEB R source files")
args <- parser$parse_args()

ctsleb_src <- resolve_ctsleb_src(args$ctsleb_src)
r_files <- list.files(ctsleb_src, pattern = "\\.[Rr]$", full.names = TRUE)
if (length(r_files) == 0) {
  stop("No CTSLEB R source files found in: ", ctsleb_src)
}
invisible(lapply(sort(r_files), source))

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

# Plink2 does not guarantee sscore row order === fam order; align on IID.
fam <- fread(paste0(args$target_test_plink, ".fam"), header = FALSE)
fam_iid <- as.character(fam[[2]])

align_fam_order <- function(x) {
  iid_hits <- vapply(seq_len(min(4, ncol(x))),
                     function(j) sum(as.character(x[[j]]) %in% fam_iid),
                     integer(1))
  iid_col <- which.max(iid_hits)
  ord <- match(fam_iid, as.character(x[[iid_col]]))
  if (anyNA(ord)) {
    stop("CTSLEB: could not align PRS matrix to study .fam (missing IIDs in sscore).")
  }
  x[ord, ]
}

prs_mat <- align_fam_order(prs_mat)

read_y <- function(f) {
  yy <- fread(f)
  if (ncol(yy) >= 2) {
    as.numeric(yy[[2]])
  } else {
    as.numeric(yy[[1]])
  }
}
y_tune <- read_y(args$tuning_pheno)
y_valid <- read_y(args$validation_pheno)
n_tune <- length(y_tune)

if (nrow(prs_mat) < (n_tune + length(y_valid))) {
  stop("CTSLEB wrapper: PRS matrix rows (", nrow(prs_mat),
       ") are fewer than tuning + validation samples (", n_tune + length(y_valid), ").")
}

prs_tune <- prs_mat[seq_len(n_tune), ]
prs_valid <- prs_mat[seq.int(n_tune + 1, n_tune + length(y_valid)), ]

# Tuning R^2 across all 2D CT PRS columns; pick the best SNP set.
prs_cols <- seq.int(3, ncol(prs_tune))
prs_r2 <- vapply(prs_cols, function(i) {
  summary(lm(y_tune ~ prs_tune[[i]]))$r.squared
}, numeric(1))
max_ind <- which.max(replace(prs_r2, !is.finite(prs_r2), -Inf))
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
prs_mat_eb <- align_fam_order(prs_mat_eb)

prs_tune_eb <- prs_mat_eb[seq_len(n_tune), ]
prs_valid_eb <- prs_mat_eb[seq.int(n_tune + 1, n_tune + length(y_valid)), ]

# Step 3a: clean the EB PRS matrices for the super-learner.
cleaned <- PRS_Clean(
  Tune_PRS = prs_tune_eb,
  Tune_Y = data.frame(y = y_tune),
  Validation_PRS = prs_valid_eb
)
clean_tune <- cleaned$Cleaned_Tune_PRS
clean_valid <- cleaned$Cleaned_Validation_PRS

r2_file <- file.path(args$results_dir, paste0(args$out_prefix, "_sl_validation_r2.txt"))
write_metric <- function(r2, note) {
  lines <- c(
    paste0("# CT-SLEB super-learner validation R2 (out_prefix=", args$out_prefix, ")"),
    "validation_r2" = paste(r2)
  )
  if (!is.null(note)) {
    lines <- c(lines, paste0("# note: ", note))
  }
  writeLines(lines, r2_file)
}

n_prs_full <- ncol(clean_tune) - 2
if (n_prs_full < 1) {
  write_metric("NA", "All PRS columns dropped during PRS_Clean; no predictors remain for the super-learner.")
  cat("CTSLEB: no super-learner predictors left after cleaning.\n")
  q(save = "no", status = 0)
}

X_tune_sl <- as.data.frame(clean_tune[, -(1:2)])
X_valid_sl <- as.data.frame(clean_valid[, -(1:2)])

# Drop near-constant columns (from the tuning partition) so glmnet/SL are stable.
drops <- caret::nearZeroVar(X_tune_sl)
if (length(drops) > 0) {
  X_tune_sl <- X_tune_sl[, -drops, drop = FALSE]
  X_valid_sl <- X_valid_sl[, -drops, drop = FALSE]
}
if (ncol(X_tune_sl) == 0) {
  write_metric("NA", "All PRS columns were near-constant after cleaning; no predictors remain.")
  cat("CTSLEB: all cleaned PRS columns near-constant; skipping super-learner.\n")
  q(save = "no", status = 0)
}

# Step 3b: super-learner on the tuning partition, evaluate on validation.
sl <- SuperLearner(
  Y = y_tune,
  X = X_tune_sl,
  family = gaussian(),
  SL.library = c("SL.glmnet", "SL.ridge")
)
y_pred_valid <- predict(sl, X_valid_sl, onlySL = TRUE)$pred
r2_valid <- summary(lm(y_valid ~ y_pred_valid))$r.squared
write_metric(r2_valid, NULL)

# Step 4: per-SNP final coefficients (SuperLearner-weighted SNP effects).
y_pred_tune <- predict(sl, X_tune_sl, onlySL = TRUE)$pred
final_betas <- ExtractFinalBetas(
  Tune_PRS = X_tune_sl,
  Predicted_Tune_Y = y_pred_tune,
  prs_mat_eb = prs_mat_eb,
  unique_infor_post = unique_infor_post,
  pthres = pthres
)
fwrite(as.data.table(final_betas),
       file.path(args$results_dir, paste0(args$out_prefix, "_final_coefficients.txt")),
       sep = "\t")

# Intermediates (also informative / debuggable).
fwrite(as.data.table(prs_mat),
       file.path(args$results_dir, paste0(args$out_prefix, "_dimct_prs.tsv")), sep = "\t")
fwrite(as.data.table(prs_mat_eb),
       file.path(args$results_dir, paste0(args$out_prefix, "_eb_prs.tsv")), sep = "\t")
fwrite(as.data.table(clean_tune),
       file.path(args$results_dir, paste0(args$out_prefix, "_tune_clean.tsv")), sep = "\t")
fwrite(as.data.table(clean_valid),
       file.path(args$results_dir, paste0(args$out_prefix, "_validation_clean.tsv")), sep = "\t")
fwrite(data.table(best_snps = best_snps),
       file.path(args$results_dir, paste0(args$out_prefix, "_best_snps.tsv")), sep = "\t")

cat("CT-SLEB wrapper completed.\n")
cat("Best SNP set:", best_snps, "\n")
cat("Super-learner validation R2:", r2_valid, "\n")