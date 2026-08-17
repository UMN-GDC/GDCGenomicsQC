#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(jsonlite)
  library(splines)
})

parser <- ArgumentParser(description = "Simulate a longitudinal HM3 phenotype with a genetic spline trajectory")
parser$add_argument("--fam", required = TRUE)
parser$add_argument("--score", required = TRUE)
parser$add_argument("--out", required = TRUE)
parser$add_argument("--basis-out", required = TRUE)
parser$add_argument("--visits", required = TRUE)
parser$add_argument("--basis-df", type = "integer", default = 4)
parser$add_argument("--heritability", type = "double", default = 0.4)
parser$add_argument("--residual-rho", type = "double", default = 0.5)
parser$add_argument("--mean-intercept", type = "double", default = 0)
parser$add_argument("--mean-slope", type = "double", default = 0)
parser$add_argument("--seed", type = "integer", default = 42)
args <- parser$parse_args()

if (args$heritability <= 0 || args$heritability >= 1) {
  stop("--heritability must be between 0 and 1")
}
if (args$residual_rho < 0 || args$residual_rho > 1) {
  stop("--residual-rho must be between 0 and 1")
}

set.seed(args$seed)
visits <- as.numeric(jsonlite::fromJSON(args$visits))
if (length(visits) <= args$basis_df) {
  stop("The number of visits must be greater than --basis-df")
}

basis <- splines::ns(visits, df = args$basis_df, intercept = TRUE)
colnames(basis) <- paste0("B", seq_len(ncol(basis)))
basis_dt <- data.table(time = visits, basis)
fwrite(basis_dt, args$basis_out, sep = "\t")

fam <- fread(args$fam, header = FALSE)
setnames(fam, c("FID", "IID", "PID", "MID", "SEX", "PHENO"))
fam <- fam[, .(FID = as.character(FID), IID = as.character(IID))]

scores <- fread(args$score)
iid_col <- intersect(c("IID", "#IID"), names(scores))[1]
fid_col <- intersect(c("FID", "#FID"), names(scores))[1]
score_cols <- grep("^BETA_B[0-9]+_SUM$", names(scores), value = TRUE)
if (is.na(iid_col) || length(score_cols) != args$basis_df) {
  stop("Expected IID and exactly ", args$basis_df, " *_SUM score columns in ", args$score)
}
if (is.na(fid_col)) {
  scores[, FID := get(iid_col)]
  fid_col <- "FID"
}

scores <- scores[, c(fid_col, iid_col, score_cols), with = FALSE]
setnames(scores, c("FID", "IID", paste0("G", seq_along(score_cols))))
scores[, FID := as.character(FID)]
scores[, IID := as.character(IID)]
dt <- merge(fam, scores, by = c("FID", "IID"), all = FALSE)

gmat <- as.matrix(dt[, paste0("G", seq_len(args$basis_df)), with = FALSE])
gmat <- scale(gmat)
shared_sd <- sqrt((1 - args$heritability) * args$residual_rho)
visit_sd <- sqrt((1 - args$heritability) * (1 - args$residual_rho))
shared_error <- rnorm(nrow(dt), sd = shared_sd)

long <- rbindlist(lapply(seq_along(visits), function(i) {
  raw_genetic <- as.numeric(gmat %*% basis[i, ])
  genetic <- as.numeric(scale(raw_genetic))
  mean_value <- args$mean_intercept + args$mean_slope * visits[i]
  phenotype <- mean_value + sqrt(args$heritability) * genetic + shared_error +
    rnorm(nrow(dt), sd = visit_sd)

  data.table(
    FID = dt$FID,
    IID = dt$IID,
    visit = i,
    time = visits[i],
    pheno = phenotype,
    genetic_trajectory = genetic
  )
}))

fwrite(long, args$out, sep = "\t")
cat("Wrote spline basis:", args$basis_out, "\n")
cat("Wrote longitudinal spline phenotype:", args$out, "\n")
