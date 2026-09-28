#!/usr/bin/env Rscript
# generate_test_data.R — Create test PLINK and GWAS data for pipeline testing
#
# Usage:
#   Rscript generate_test_data.R <out_dir> [N] [M] [--seed=<n>]
#
# Generates:
#   PLINK trio: study_sample.{bed,bim,fam} + .bk + .rds (bigSNP)
#   GWAS files: gwas_pheno{A,B,C}.txt (raw format, 3 columns: SNP BETA P)
#   GWAS files: aligned_gwas_{phenoX}.txt (pre-aligned format with full header)
#   Combined: combined_phenotypes.txt (FID IID phenoA phenoB phenoC)
#   Combined: study_sample_pheno.txt (FID IID phenoX — one per phenotype)
#   Configs: test_config.txt, test_config_skip_ss.txt

args <- commandArgs(trailingOnly = TRUE)

out_dir <- args[1]
if (is.na(out_dir)) stop("Usage: generate_test_data.R <out_dir> [N] [M] [--seed=<n>]")

N <- as.integer(args[2] %||% "100")
M <- as.integer(args[3] %||% "200")

seed <- 42
if (any(grepl("^--seed=", args))) {
  seed <- as.integer(sub("--seed=", "", grep("^--seed=", args, value = TRUE)[1]))
}
set.seed(seed)

`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- 1. Simulate SNP info ---
chromosomes <- sample(rep(1:22, length.out = M))
positions <- sort(sample(1:1e6, M))
snp_ids <- sprintf("rs%d", seq_len(M))
allele1 <- sample(c("A", "C", "G", "T"), M, replace = TRUE)
allele2 <- sapply(allele1, function(a1) {
  rest <- setdiff(c("A", "C", "G", "T"), a1)
  sample(rest, 1)
})
genetic_dist <- rep(0, M)

# --- 2. Simulate genotypes (HWE-like: p ~ uniform, binomial) ---
afs <- runif(M, 0.1, 0.5)
genotypes <- matrix(NA_integer_, N, M)
for (j in seq_len(M)) {
  p <- afs[j]
  # HWE: P(0)=(1-p)^2, P(1)=2p(1-p), P(2)=p^2
  probs <- c((1-p)^2, 2*p*(1-p), p^2)
  genotypes[, j] <- sample(0:2, N, replace = TRUE, prob = probs)
}

# --- 3. FID/IID ---
fid <- rep("Test", N)
iid <- sprintf("IND%04d", seq_len(N))

# --- 4. Simulate 3 phenotypes (continuous, ~3000-5000) ---
pheno_names <- c("phenoA", "phenoB", "phenoC")
n_pheno <- length(pheno_names)

# Generate correlated phenotypes with genetic signal
genetic_effects <- rep(0, N)
for (j in 1:10) {  # first 10 SNPs are causal
  genetic_effects <- genetic_effects + 0.5 * (genotypes[, j] - 2 * afs[j])
}
genetic_effects <- scale(genetic_effects)

pheno_matrix <- matrix(NA_real_, N, n_pheno)
colnames(pheno_matrix) <- pheno_names

# Each phenotype has some genetic signal + noise, correlated
for (k in seq_len(n_pheno)) {
  noise <- rnorm(N, mean = 4000, sd = 300)
  pheno_matrix[, k] <- 4000 + 200 * genetic_effects + noise
}
pheno_matrix <- pmax(pheno_matrix, 3000)  # floor at 3000

# --- 5. Write FAM file ---
fam <- data.frame(
  family.ID = fid,
  sample.ID = iid,
  paternal.ID = rep(0, N),
  maternal.ID = rep(0, N),
  sex = rep(1, N),
  affection = rep(-9, N)
)
utils::write.table(fam, file.path(out_dir, "study_sample.fam"),
  row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# --- 6. Write BIM file ---
bim <- data.frame(
  chr = chromosomes,
  rsid = snp_ids,
  gd = genetic_dist,
  bp = positions,
  a1 = allele1,
  a2 = allele2
)
utils::write.table(bim, file.path(out_dir, "study_sample.bim"),
  row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# --- 7. Write MAP file (for C+T / PRSice2 text-mode) ---
map <- data.frame(chr = chromosomes, rsid = snp_ids, gd = genetic_dist, bp = positions)
utils::write.table(map, file.path(out_dir, "study_sample.map"),
  row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# --- 8. Write PED file (for C+T / PRSice2 text mode) ---
# Format: FID IID father mother sex pheno + space-sep alleles (A1 A2 pairs)
ped_df <- data.frame(
  FID = fid, IID = iid, father = 0, mother = 0, sex = 1, pheno = -9,
  stringsAsFactors = FALSE
)
gt_alleles <- apply(genotypes, 1, function(gt_row) {
  paste(sapply(gt_row, function(g) {
    switch(as.character(g),
      "0" = paste(allele2, allele2),
      "1" = paste(allele1, allele2),
      "2" = paste(allele1, allele1),
      "0 0")
  }), collapse = " ")
})
ped_df$alleles <- gt_alleles
writeLines(apply(ped_df, 1, paste, collapse = " "),
  file.path(out_dir, "study_sample.ped"))

message("Wrote: study_sample.{fam,bim,ped}")

# --- 9. Create bigSNP .bk + .rds via snp_fake() ---
message("Creating bigSNP via snp_fake()...")
suppressPackageStartupMessages(library(bigsnpr))

fake <- snp_fake(n = N, m = M)
fake$fam <- data.frame(
  family.ID = fid, sample.ID = iid,
  paternal.ID = rep(0, N), maternal.ID = rep(0, N),
  sex = rep(1, N), affection = rep(NA_real_, N)
)
fake$map <- data.frame(
  chromosome = chromosomes, marker.ID = snp_ids,
  genetic.dist = genetic_dist, physical.pos = positions,
  allele1 = allele1, allele2 = allele2,
  stringsAsFactors = FALSE
)

# Fill genotypes
fake$genotypes[, ] <- genotypes
snp_save(fake)

# Copy .bk and .rds to out_dir
out_bk <- file.path(out_dir, "study_sample.bk")
out_rds <- file.path(out_dir, "study_sample.rds")
file.copy(fake$genotypes$bk, out_bk, overwrite = TRUE)
file.copy(fake$genotypes$rds, out_rds, overwrite = TRUE)

# Verify
verify <- snp_attach(out_rds)
cat(sprintf("  bigSNP: %d x %d, G[1,1]=%d\n",
  nrow(verify$genotypes), ncol(verify$genotypes), verify$genotypes[1, 1]))

rm(fake, verify)
invisible(gc())

# --- 10. Write GWAS summary stats (raw format: SNP BETA P) ---
set.seed(seed + 1)
gwas_dir <- file.path(out_dir, "gwas")
dir.create(gwas_dir, showWarnings = FALSE)

# Use only independent SNPs for GWAS signal (prune to 50 causal-like SNPs)
prune_idx <- sort(sample(M, min(50, M)))  # pick 50 representative SNPs
nw <- M  # keep all SNPs in output but only prune_idx have strong signal

for (k in seq_len(n_pheno)) {
  pn <- pheno_names[k]
  pheno <- pheno_matrix[, k]

  # Compute crude betas by regressing phenotype on each SNP
  beta <- se <- pval <- rep(NA_real_, M)
  for (j in seq_len(M)) {
    g <- genotypes[, j]
    keep <- !is.na(g)
    if (sum(keep) < 10) next
    y <- pheno[keep]
    x <- g[keep]
    fit <- tryCatch(lm(y ~ x), error = function(e) NULL)
    if (is.null(fit)) next
    coefs <- summary(fit)$coefficients
    if (nrow(coefs) < 2) next
    beta[j] <- coefs[2, "Estimate"]
    se[j] <- coefs[2, "Std. Error"]
    pval[j] <- coefs[2, "Pr(>|t|)"]
  }
  # Replace NAs with small-signal noise
  beta[is.na(beta)] <- rnorm(sum(is.na(beta)), 0, 0.5)
  se[is.na(se)] <- abs(rnorm(sum(is.na(se)), 40, 10))
  pval[is.na(pval)] <- runif(sum(is.na(pval)))

  # Raw format: 3 columns, no header (for prepare_sumstats.R)
  raw_file <- file.path(gwas_dir, sprintf("gwas_%s.txt", pn))
  raw_df <- data.frame(SNP = snp_ids, BETA = beta, P = pval)
  utils::write.table(raw_df, raw_file, row.names = FALSE, col.names = TRUE,
    quote = FALSE, sep = "\t")
  message(sprintf("  Wrote: gwas/%s (raw, N=%d)", basename(raw_file), M))
}

# --- 11. Write pre-aligned GWAS for skip_ss mode ---
for (k in seq_len(n_pheno)) {
  pn <- pheno_names[k]
  raw_file <- file.path(gwas_dir, sprintf("gwas_%s.txt", pn))
  raw <- data.table::fread(raw_file)

  # Merge with BIM to get CHR, BP, A1, A2, add beta_se, n_eff
  raw <- merge(raw, bim, by.x = "SNP", by.y = "rsid", all.x = TRUE, sort = FALSE)
  raw$beta_se <- abs(raw$BETA) / max(1, qnorm(raw$P / 2 + 1e-16))  # approximate SE
  raw$beta_se[is.infinite(raw$beta_se) | is.na(raw$beta_se)] <- 50
  raw$n_eff <- N

  aligned <- raw[, c("SNP", "chr", "bp", "a1", "a2", "BETA", "beta_se", "P", "n_eff")]
  names(aligned) <- c("SNP", "CHR", "BP", "A1", "A2", "beta", "beta_se", "P", "n_eff")

  aligned_file <- file.path(gwas_dir, sprintf("aligned_gwas_%s.txt", pn))
  utils::write.table(aligned, aligned_file, row.names = FALSE, col.names = TRUE,
    quote = FALSE, sep = "\t")
  message(sprintf("  Wrote: gwas/%s (aligned)", basename(aligned_file)))
}

# --- 12. Write combined phenotype file ---
pheno_df <- data.frame(FID = fid, IID = iid, pheno_matrix)
utils::write.table(pheno_df, file.path(out_dir, "combined_phenotypes.txt"),
  row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
message("  Wrote: combined_phenotypes.txt")

# Also write per-phenotype files
for (k in seq_len(n_pheno)) {
  pn <- pheno_names[k]
  pp <- file.path(gwas_dir, sprintf("study_sample_pheno_%s.txt", pn))
  pp_df <- data.frame(FID = fid, IID = iid, pheno_matrix[, k])
  utils::write.table(pp_df, pp, row.names = FALSE, col.names = TRUE,
    quote = FALSE, sep = "\t")
}
message("  Wrote: gwas/study_sample_pheno_*.txt per phenotype")

# --- 13. Write config files ---
config <- c(
  sprintf('output_path="%s"', out_dir),
  "# Multiple summary stats (comma-separated)",
  sprintf('summary_stats_files="%s"', paste(
    file.path(gwas_dir, sprintf("aligned_gwas_%s.txt", pheno_names)),
    collapse = ","
  )),
  sprintf('phenotype_info_file="%s/combined_phenotypes.txt"', out_dir),
  sprintf('anc_PREFIX="%s/study_sample"', out_dir),
  sprintf('anc_bed="%s/study_sample"', out_dir),
  'ref_1000G_prefix=""',
  'ncores=2',
  'run_ct=1',
  'run_ldpred2=1',
  'run_lassosum2=1',
  'run_prsice2=1',
  'ld_window=1000',
  'ld_r2=0.1',
  'pval_thresholds=5e-8,5e-7,5e-6,5e-5,5e-4,5e-3,5e-2,0.1,0.2,0.5,1.0',
  'prslice_r2=0.1',
  'prslice_window=1000',
  'prslice_pval=0.05,0.1,0.5,1.0',
  'ldpred2_model=grid',
  'ldpred2_pval=0.05,0.1,0.5,1.0',
  'lassosum2_pval=0.05,0.1,0.5,1.0',
  'snp_ids_col=SNP',
  'beta_col=beta',
  'pval_col=P',
  'n_total=100',
  'plink_allowed_chr=1-22',
  'pca_include=T',
  'grm_include=F',
  'split_grm=F',
  paste0('pcs=', paste0("PC", 1:10, collapse = ","))
)

writeLines(config, file.path(out_dir, "test_config.txt"))
message("  Wrote: test_config.txt")

# Full config but skip_ss_generation=1 (for faster iteration)
config_ss <- sub("summary_stats_files=.*",
  paste(sprintf('summary_stats_files="%s"', paste(
    file.path(gwas_dir, sprintf("aligned_gwas_%s.txt", pheno_names)),
    collapse = ","
  ))),
  config)
config_ss <- sub("^run_ct=.*", "run_ct=0", config_ss)
config_ss <- sub("^run_ldpred2=.*", "run_ldpred2=1", config_ss)
config_ss <- sub("^run_lassosum2=.*", "run_lassosum2=0", config_ss)
config_ss <- sub("^run_prsice2=.*", "run_prsice2=0", config_ss)
writeLines(config_ss, file.path(out_dir, "test_config_ldpred2.txt"))

# Single-phenotype config for backward-compat test
config_sp <- config
config_sp[grep("^summary_stats_files=", config_sp)] <-
  sprintf('summary_stats_files="%s/gwas/aligned_gwas_phenoA.txt"', out_dir)
config_sp[grep("^phenotype_info_file=", config_sp)] <-
  sprintf('phenotype_info_file="%s/gwas/study_sample_pheno_phenoA.txt"', out_dir)
writeLines(config_sp, file.path(out_dir, "test_config_single.txt"))

message("  Wrote: test_config_single.txt (single-phenotype backward compat)")

# --- Summary ---
cat(sprintf("\n=== Test data generated at %s ===\n", out_dir))
cat(sprintf("  Individuals: %d\n", N))
cat(sprintf("  SNPs: %d\n", M))
cat(sprintf("  Phenotypes: %s\n", paste(pheno_names, collapse = ", ")))
cat(sprintf("  Pheno ranges: %.0f-%.0f\n",
  min(pheno_matrix), max(pheno_matrix)))
cat(sprintf("  PLINK: %s/study_sample.{bed,bim,fam}\n", out_dir))
cat(sprintf("  bigSNP: %s/study_sample.{rds,bk}\n", out_dir))
cat(sprintf("  GWAS: %s/gwas/aligned_gwas_{%s}.txt\n", out_dir, paste(pheno_names, collapse = ",")))
cat(sprintf("  Pheno: %s/combined_phenotypes.txt\n", out_dir))
cat(sprintf("  Configs: %s/test_config*.txt\n", out_dir))

invisible(0)
