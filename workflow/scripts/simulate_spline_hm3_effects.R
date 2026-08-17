#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(MASS)
})

parser <- ArgumentParser(description = "Simulate correlated HM3 genetic effects for a spline trajectory")
parser$add_argument("--bim", required = TRUE)
parser$add_argument("--out-prefix", required = TRUE)
parser$add_argument("--n-causal", type = "integer", default = 1000)
parser$add_argument("--basis-df", type = "integer", default = 4)
parser$add_argument("--basis-rho", type = "double", default = 0.5)
parser$add_argument("--seed", type = "integer", default = 42)
args <- parser$parse_args()

if (args$basis_df < 2) stop("--basis-df must be at least 2")
if (abs(args$basis_rho) >= 1) stop("--basis-rho must be between -1 and 1")

set.seed(args$seed)

bim <- fread(args$bim, header = FALSE)
setnames(bim, c("CHR", "SNP", "CM", "BP", "A1", "A2"))
if (args$n_causal > nrow(bim)) {
  stop("--n-causal is larger than the number of variants in the BIM file")
}

if (args$n_causal <= 0) {
  causal <- copy(bim)
} else {
  causal <- bim[sample.int(.N, args$n_causal)]
}

index <- seq_len(args$basis_df)
sigma <- outer(index, index, function(i, j) args$basis_rho^abs(i - j))
beta_names <- paste0("BETA_B", index)

nvar <- nrow(causal)
chunk_size <- 100000L
n_chunks <- ceiling(nvar / chunk_size)

effect_table <- data.table(SNP = causal$SNP, A1 = causal$A1)
for (nm in beta_names) {
  effect_table[, (nm) := numeric(nvar)]
}

for (chunk_idx in seq_len(n_chunks)) {
  start <- (chunk_idx - 1L) * chunk_size + 1L
  stop <- min(chunk_idx * chunk_size, nvar)
  n_this <- stop - start + 1L

  chunk_effects <- MASS::mvrnorm(
    n = n_this,
    mu = rep(0, args$basis_df),
    Sigma = sigma
  )

  for (j in seq_len(args$basis_df)) {
    effect_table[start:stop, (beta_names[j]) := chunk_effects[, j]]
  }

  cat("Finished chunk", chunk_idx, "of", n_chunks, "\n")
}

# Keep each basis-specific polygenic score on a comparable scale.
for (nm in beta_names) {
  x <- effect_table[[nm]]
  x <- as.numeric(scale(x)) / sqrt(length(x))
  effect_table[, (nm) := x]
}



effect_file <- paste0(args$out_prefix, "_effects.tsv")
score_file <- paste0(args$out_prefix, "_basis.score")
fwrite(effect_table, effect_file, sep = "\t")
fwrite(effect_table[, c("SNP", "A1", beta_names), with = FALSE], score_file, sep = "\t")

cat("Wrote spline effect table:", effect_file, "\n")
cat("Wrote multi-column PLINK score file:", score_file, "\n")
cat("Variants with simulated effects:", nrow(causal), "\n")
