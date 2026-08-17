#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
})

parser <- ArgumentParser(description = "Derive subject-specific spline trajectory coefficients for AdjHE")
parser$add_argument("--longitudinal", required = TRUE)
parser$add_argument("--basis", required = TRUE)
parser$add_argument("--grm-id", required = TRUE)
parser$add_argument("--out", required = TRUE)
args <- parser$parse_args()

long <- fread(args$longitudinal)
basis <- fread(args$basis)
grm_ids <- fread(args$grm_id, header = FALSE)
setnames(grm_ids, c("FID", "IID"))
grm_ids[, `:=`(FID = as.character(FID), IID = as.character(IID), grm_order = .I)]

required <- c("FID", "IID", "time", "pheno")
missing <- setdiff(required, names(long))
if (length(missing)) stop("Missing longitudinal columns: ", paste(missing, collapse = ", "))

basis_cols <- grep("^B[0-9]+$", names(basis), value = TRUE)
if (!length(basis_cols)) stop("No B1...Bk columns found in basis file")

long[, `:=`(FID = as.character(FID), IID = as.character(IID))]
long <- merge(long, basis, by = "time", all.x = TRUE, sort = FALSE)
if (anyNA(long[, ..basis_cols])) stop("Some visit times did not match the spline basis file")

fit_subject <- function(dt) {
  complete <- complete.cases(dt[, c("pheno", basis_cols), with = FALSE])
  y <- dt$pheno[complete]
  x <- as.matrix(dt[complete, ..basis_cols])
  if (nrow(x) < ncol(x) || qr(x)$rank < ncol(x)) {
    return(as.list(rep(NA_real_, length(basis_cols))))
  }
  as.list(as.numeric(qr.coef(qr(x), y)))
}

coef_dt <- long[, fit_subject(.SD), by = .(FID, IID)]
setnames(coef_dt, c("FID", "IID", paste0("spline_coef_B", seq_along(basis_cols))))

coef_dt <- merge(grm_ids, coef_dt, by = c("FID", "IID"), all.x = TRUE, sort = FALSE)
setorder(coef_dt, grm_order)
coef_dt[, grm_order := NULL]
coef_cols <- setdiff(names(coef_dt), c("FID", "IID"))
coef_dt <- coef_dt[complete.cases(coef_dt[, ..coef_cols])]

if (!nrow(coef_dt)) stop("No complete spline coefficients remained after GRM-ID matching")
fwrite(coef_dt, args$out, sep = " ", col.names = FALSE)

cat("Wrote AdjHE spline coefficient phenotype:", args$out, "\n")
cat("Subjects:", nrow(coef_dt), "coefficients:", length(coef_cols), "\n")
