#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
})

parser <- ArgumentParser(description = "Prepare CTSLEB-formatted summary stats")
parser$add_argument("--sumstats", required = TRUE)
parser$add_argument("--bim", required = TRUE)
parser$add_argument("--out", required = TRUE)
args <- parser$parse_args()

ss <- fread(args$sumstats)
bim <- fread(args$bim, header = FALSE)
setnames(bim, c("CHR_bim", "SNP", "CM", "BP", "A1_bim", "A2_bim"))

x <- merge(ss, bim[, .(SNP, BP)], by = "SNP", all.x = TRUE)

required_cols <- c("CHR", "SNP", "A1", "A2", "BETA", "SE", "P")
missing_cols <- setdiff(required_cols, names(x))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

x[, rs_id := SNP]

if (!"N" %in% names(x)) {
  x[, N := NA_real_]
}

out <- x[, .(CHR, SNP, BP, A1, A2, BETA, SE, P, rs_id, N)]

if (anyNA(out$BP)) {
  stop("Some BP values are missing after BIM join.")
}

fwrite(out, args$out, sep = "\t")
cat("Wrote:", args$out, "\n")
