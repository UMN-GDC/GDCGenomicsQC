#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(MKL_NUM_THREADS = "1")

library(argparse)
library(bigsnpr)
options(bigstatsr.check.parallel = FALSE)

parser <- ArgumentParser(description="Pre-compute LD matrices from genotype data for reuse across PRS methods")
parser$add_argument("--anc_bed", type="character", required=TRUE,
                    help="Path to PLINK .bed file")
parser$add_argument("--out", type="character", required=TRUE,
                    help="Output directory (e.g. gwas/ld_matrix/)")
parser$add_argument("--ncores", type="integer", default=1,
                    help="Number of CPU cores")
args <- parser$parse_args()

dir.create(args$out, recursive = TRUE, showWarnings = FALSE)

# Load genotype data
rds_path <- sub("\\.bed$", ".rds", args$anc_bed)
if (!file.exists(rds_path)) {
  message("Converting .bed to .rds...")
  snp_readBed(args$anc_bed)
}
obj.bigSNP <- snp_attach(rds_path)
G <- obj.bigSNP$genotypes
POS2 <- obj.bigSNP$map$genetic.dist
n_total <- nrow(G)

message(sprintf("Loaded: %d samples, %d variants", n_total, ncol(G)))

CHRS <- 1:22
saved_chrs <- integer()
total_snps <- 0

for (chr in CHRS) {
  ind.chr <- which(obj.bigSNP$map$chromosome == chr)
  n_chr <- length(ind.chr)
  if (n_chr < 2) next

  message(sprintf("Chr%d: %d SNPs ... computing", chr, n_chr))

  corr0 <- snp_cor(G, ind.col = ind.chr, size = 3 / 1000,
                   infos.pos = POS2[ind.chr], ncores = args$ncores)

  # Repair any NA/NaN correlations (monomorphic or fully-missing SNPs)
  if (any(is.na(corr0))) {
    na_count <- sum(is.na(corr0))
    message(sprintf("  Repaired %d NA correlations", na_count))
    corr0[is.na(corr0)] <- 0
  }
  # Add small ridge to ensure positive definiteness for downstream Gibbs samplers
  corr0 <- corr0 + diag(ncol(corr0)) * 1e-5

  saveRDS(corr0, file.path(args$out, sprintf("chr%d_corr.rds", chr)))
  message(sprintf("  Saved chr%d_corr.rds", chr))

  # Free memory — correlation matrices can exceed 500MB for large chromosomes
  rm(corr0, ind.chr)
  gc()

  saved_chrs <- c(saved_chrs, chr)
  total_snps <- total_snps + n_chr
}

# Save full map (chr, rsid, pos, a1, a0) matching the LD matrix order
# obj.bigSNP$map columns: chromosome, rsid, genetic.dist, physical.pos, a1, a0
map_full <- obj.bigSNP$map[, c(1, 2, 4, 5, 6)]
colnames(map_full) <- c("chr", "rsid", "pos", "a1", "a0")
map_full <- map_full[map_full$chr %in% saved_chrs, ]
saveRDS(map_full, file.path(args$out, "map.rds"))

# Save mapping from map row -> original G column index
g_idx <- which(obj.bigSNP$map$chromosome %in% saved_chrs)
saveRDS(g_idx, file.path(args$out, "g_idx.rds"))

message(sprintf("\nDone! %d total SNPs saved to %s", total_snps, args$out))
