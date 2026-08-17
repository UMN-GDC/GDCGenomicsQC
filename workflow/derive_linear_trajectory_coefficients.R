#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(jsonlite)
})

parser <- ArgumentParser(description = "Derive subject-specific linear trajectory coefficients from longitudinal simulation data")
parser$add_argument("--longitudinal", required = TRUE)
parser$add_argument("--grm-id", required = TRUE)
parser$add_argument("--out", required = TRUE)
parser$add_argument("--manifest", required = TRUE)
parser$add_argument("--center-time", type = "double", default = 0)
args <- parser$parse_args()

long <- fread(args$longitudinal)
grm_ids <- fread(args$grm_id, header = FALSE)
setnames(grm_ids, c("FID", "IID"))
grm_ids[, `:=`(FID = as.character(FID), IID = as.character(IID), grm_order = .I)]

required <- c("FID", "IID", "time", "pheno")
missing <- setdiff(required, names(long))
if (length(missing)) stop("Missing longitudinal columns: ", paste(missing, collapse = ", "))

long[, `:=`(FID = as.character(FID), IID = as.character(IID))]
setorder(long, FID, IID, time)

fit_subject <- function(dt, center_time) {
  x <- dt$time - center_time
  y <- dt$pheno
  if (length(unique(x)) < 2) {
    return(list(intercept_t0 = NA_real_, annual_slope = NA_real_, n_visits = length(y)))
  }
  fit <- lm(y ~ x)
  list(
    intercept_t0 = unname(coef(fit)[["(Intercept)"]]),
    annual_slope = unname(coef(fit)[["x"]]),
    n_visits = length(y)
  )
}

coef_dt <- long[, fit_subject(.SD, args$center_time), by = .(FID, IID)]
coef_dt <- merge(grm_ids, coef_dt, by = c("FID", "IID"), all.x = TRUE, sort = FALSE)
setorder(coef_dt, grm_order)
coef_dt[, grm_order := NULL]
coef_dt <- coef_dt[complete.cases(coef_dt[, .(intercept_t0, annual_slope)])]

if (!nrow(coef_dt)) stop("No complete linear trajectory coefficients remained after GRM-ID matching")
fwrite(coef_dt[, .(FID, IID, intercept_t0, annual_slope)], args$out, sep = " ", col.names = FALSE)

write_json(
  list(
    phenotype_names = list("intercept_t0", "annual_slope"),
    center_time = args$center_time,
    n_final = nrow(coef_dt)
  ),
  args$manifest,
  pretty = TRUE,
  auto_unbox = TRUE
)

cat("Wrote linear trajectory phenotype:", args$out, "\n")
cat("Subjects:", nrow(coef_dt), "\n")
