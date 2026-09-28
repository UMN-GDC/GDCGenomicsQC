#!/usr/bin/env Rscript

# Control BLAS threads to avoid nested parallelism with bigstatsr's multicore
Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(MKL_NUM_THREADS = "1")

library(argparse)
library(bigsnpr)
library(ggplot2)
library(dplyr)
library(bigreadr)

options(bigstatsr.check.parallel = FALSE)
options(bigstatsr.check.args = FALSE)

# 1. SET UP ARGUMENT PARSER
parser <- ArgumentParser(description='lassosum2 Pipeline for Polygenic Risk Scores')
parser$add_argument("--anc_bed", type="character", help="Path to the plink .bed file")
parser$add_argument("--rds", type="character", help="Path to the .rds file")
parser$add_argument("--ss", type="character", required=TRUE, help="Path to GWAS summary statistics")
parser$add_argument("--bim", type="character", required=TRUE, help="Path to the .bim file")
parser$add_argument("--beta_se", type="character", help="Optional: separate table with beta_se/SE")
parser$add_argument("--afreq", type="character", help="Path to PLINK2 .afreq file")
parser$add_argument("--pheno", type="character", help="Path to phenotype file (FID, IID, phenotype).")
parser$add_argument("--n_val", type="integer", default=49, help="Number of samples for validation")
parser$add_argument("--seed", type="integer", default=1, help="Seed for validation/test split")
parser$add_argument("--out", type="character", default="lassosum_out", help="Prefix for output files")
parser$add_argument("--ncores", type="integer", default=1, help="Number of CPU cores (default: 1)")
parser$add_argument("--ld-cache-dir", type="character", help="Directory to cache/reuse per-chromosome LD matrices")
parser$add_argument("--ld-matrix-dir", type="character", help="Directory with pre-computed LD matrix")

args <- parser$parse_args()

tryCatch({

# --- 2. DATA LOADING & SAMPLE MATCHING ---
if (!is.null(args$anc_bed)) {
  rds_path <- if (!is.null(args$rds)) args$rds else sub("\\.bed$", ".rds", args$anc_bed)
  if (!file.exists(rds_path)) {
    snp_readBed(args$anc_bed)
  }
  args$rds <- rds_path
}

obj.bigSNP <- snp_attach(args$rds)
G          <- obj.bigSNP$genotypes
y          <- obj.bigSNP$fam$affection
NCORES     <- args$ncores

# Map for target study sample
target_map <- setNames(obj.bigSNP$map[-3], c("chr", "rsid", "pos", "a1", "a0"))

message("Imputing missing genotype calls in study sample...")
snp_fastImputeSimple(G, method = "mean0", ncores = NCORES)

if (!is.null(args$pheno)) {
  message("Reading phenotype from: ", args$pheno)
  pheno <- fread2(args$pheno)
  
  fam <- obj.bigSNP$fam
  fam_id <- paste(fam$family.ID, fam$sample.ID, sep = "_")
  pheno_id <- paste(pheno[[1]], pheno[[2]], sep = "_")
  
  m <- match(fam_id, pheno_id)
  matched_count <- sum(!is.na(m))
  
  message("Successfully matched ", matched_count, " / ", length(fam_id), " individuals to phenotype data.")
  
  if (matched_count == 0) {
    stop("Error: Zero individuals matched between genotype .fam file and phenotype file.")
  }
  
  y <- as.numeric(pheno[[3]][m])
}

if (!is.null(y)) {
  y[y %in% c(-9, -999)] <- NA
}

# --- 3. SUMMARY STATS PREPARATION ---
qassoc <- bigreadr::fread2(args$ss) 

sumstats_mapped <- qassoc %>%
  rename(
    rsid  = any_of(c("SNP", "rsid", "ID")),
    pos   = any_of(c("BP", "pos", "position")),
    chr   = any_of(c("CHR", "chr", "chromosome")),
    a1    = any_of(c("A1", "a1", "allele1")),
    beta  = any_of(c("BETA", "beta", "eff")),
    n_eff = any_of(c("n_eff", "NMISS", "N"))
  )

if ("TEST" %in% colnames(sumstats_mapped)) {
  sumstats_mapped <- sumstats_mapped %>% filter(TEST == "ADD")
}

bim.file <- bigreadr::fread2(args$bim, select = c(1, 2, 4, 5, 6))
colnames(bim.file) <- c("chr", "bim.rsid", "pos", "bim.a1", "bim.a0")

if (!is.null(args$beta_se)) {
  beta.se.tab <- bigreadr::fread2(args$beta_se) %>% 
    rename(rsid = any_of(c("SNP", "rsid")), 
           beta_se = any_of(c("SE", "beta_se"))) %>%
    select(rsid, beta_se)
  
  sumstats <- sumstats_mapped %>% inner_join(beta.se.tab, by = "rsid")
} else {
  sumstats <- sumstats_mapped %>% rename(beta_se = any_of(c("SE", "beta_se")))
}

sumstats <- sumstats %>%
  inner_join(bim.file, by = c("chr", "pos")) %>%
  mutate(a0 = ifelse(a1 == bim.a1, bim.a0, bim.a1)) %>%
  mutate(rsid = ifelse(is.na(rsid) | rsid == "", bim.rsid, rsid)) %>%
  distinct(rsid, .keep_all = TRUE) %>%
  select(rsid, chr, pos, a1, a0, beta, beta_se, n_eff)

sumstats_clean <- sumstats[, c("rsid", "chr", "pos", "a1", "a0", "beta", "beta_se", "n_eff")]

# --- 4. MATCH SUMSTATS TO LD REFERENCE MATRIX ---
if (!is.null(args$ld_matrix_dir)) {
  message("Aligning sumstats with pre-computed LD matrix map...")
  ld_map <- readRDS(file.path(args$ld_matrix_dir, "map.rds"))
  
  m_rsid <- tryCatch(snp_match(sumstats_clean, ld_map, join_by_pos = FALSE), error = function(e) NULL)
  m_pos  <- tryCatch(snp_match(sumstats_clean, ld_map, join_by_pos = TRUE), error = function(e) NULL)
  
  n_rsid <- if (!is.null(m_rsid)) nrow(m_rsid) else 0
  n_pos  <- if (!is.null(m_pos))  nrow(m_pos)  else 0
  
  if (n_rsid >= n_pos && n_rsid > 0) {
    df_beta <- m_rsid
    message("Matched ", n_rsid, " variants with LD reference panel by RSID.")
  } else if (n_pos > 0) {
    df_beta <- m_pos
    message("Matched ", n_pos, " variants with LD reference panel by Position.")
  } else {
    stop("Error: Failed to match summary statistics with LD reference panel.")
  }
  
  if (!is.null(args$afreq)) {
    afreq <- fread2(args$afreq)
    m <- match(df_beta$rsid, afreq$ID)
    maf <- pmin(afreq$ALT_FREQS[m], 1 - afreq$ALT_FREQS[m])
  } else if ("af" %in% names(ld_map)) {
    ld_af <- ld_map$af[df_beta$`_NUM_ID_`]
    maf <- pmin(ld_af, 1 - ld_af)
  } else {
    maf <- rep(NA_real_, nrow(df_beta))
  }
} else {
  ld_map <- target_map
  df_beta <- tryCatch({
    snp_match(sumstats_clean, target_map, join_by_pos = FALSE)
  }, error = function(e) {
    snp_match(sumstats_clean, target_map, join_by_pos = TRUE)
  })
  
  if (!is.null(args$afreq)) {
    afreq <- fread2(args$afreq)
    m <- match(df_beta$rsid, afreq$ID)
    maf <- pmin(afreq$ALT_FREQS[m], 1 - afreq$ALT_FREQS[m])
  } else {
    maf <- snp_MAF(G, ind.col = df_beta$`_NUM_ID_`, ncores = NCORES)
  }
}

maf_thr <- 1 / sqrt(nrow(G))
keep_maf <- is.na(maf) | (maf > maf_thr)
df_beta  <- df_beta[keep_maf, ]

CHRS <- 1:22
df_beta <- df_beta[df_beta$chr %in% CHRS, ]
message("   Autosomes 1-22 retained for modeling: ", nrow(df_beta))

if (nrow(df_beta) == 0) {
  stop("No variants remain after MAF + autosome filtering.")
}

# --- 5. BUILD LD CORRELATION MATRIX (corr) ---
if (!is.null(args$ld_cache_dir)) {
  dir.create(args$ld_cache_dir, showWarnings = FALSE, recursive = TRUE)
}

tmp <- tempfile(tmpdir = "temp_ld_lassosum")
dir.create("temp_ld_lassosum", showWarnings = FALSE)
corr <- NULL
df_beta_list <- list()

if (!is.null(args$ld_matrix_dir)) {
  message("Loading pre-computed LD matrix from: ", args$ld_matrix_dir)

  for (chr in CHRS) {
    ind.chr <- which(df_beta$chr == chr)
    if (length(ind.chr) < 2) next

    ld_map_chr_idx <- which(ld_map$chr == chr)
    local_idx <- match(df_beta$`_NUM_ID_`[ind.chr], ld_map_chr_idx)
    
    bad <- which(is.na(local_idx))
    if (length(bad) > 0) {
      local_idx <- local_idx[-bad]
      ind.chr <- ind.chr[-bad]
    }
    if (length(local_idx) < 2) next

    possible_files <- c(
      file.path(args$ld_matrix_dir, sprintf("chr%d_corr.rds", chr)),
      file.path(args$ld_matrix_dir, sprintf("chr%d.rds", chr)),
      file.path(args$ld_matrix_dir, sprintf("LD_chr%d.rds", chr))
    )
    corr_file <- possible_files[file.exists(possible_files)][1]

    if (is.na(corr_file) || !file.exists(corr_file)) {
      message("Warning: LD matrix file not found for chromosome ", chr, ". Skipping.")
      next
    }

    corr0_full <- readRDS(corr_file)
    corr0 <- corr0_full[local_idx, local_idx]

    if (any(is.na(corr0))) corr0[is.na(corr0)] <- 0
    corr0 <- corr0 + diag(ncol(corr0)) * 1e-5

    if (is.null(corr)) {
      corr <- as_SFBM(corr0, tmp, compact = TRUE)
    } else {
      corr$add_columns(corr0, nrow(corr))
    }
    df_beta_list[[length(df_beta_list) + 1]] <- df_beta[ind.chr, ]
  }
  df_beta <- do.call(rbind, df_beta_list)
  df_beta$`_NUM_ID_` <- seq_len(nrow(df_beta))
} else {
  keep_idx <- logical(nrow(df_beta))
  POS2 <- obj.bigSNP$map$genetic.dist
  for (chr in CHRS) {
    ind.chr <- which(df_beta$chr == chr)
    ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
    if (length(ind.chr2) < 2) next

    keep_idx[ind.chr] <- TRUE
    corr0 <- snp_cor(G, ind.col = ind.chr2, size = 3/1000, infos.pos = POS2[ind.chr2], ncores = NCORES)

    if (any(is.na(corr0))) corr0[is.na(corr0)] <- 0
    corr0 <- corr0 + diag(ncol(corr0)) * 1e-5

    if (is.null(corr)) {
      corr <- as_SFBM(corr0, tmp, compact = TRUE)
    } else {
      corr$add_columns(corr0, nrow(corr))
    }
    df_beta_list[[length(df_beta_list) + 1]] <- df_beta[ind.chr, ]
  }
  df_beta <- do.call(rbind, df_beta_list)
  df_beta$`_NUM_ID_` <- seq_len(nrow(df_beta))
}

message("Variants in LD matrix: ", nrow(df_beta))

if (nrow(df_beta) != ncol(corr)) {
  stop("Mismatch between df_beta rows (", nrow(df_beta), ") and corr cols (", ncol(corr), ") after LD construction.")
}

# --- 6. LASSOSUM2 MODELING ---
beta_lassosum2 <- snp_lassosum2(corr, df_beta, ncores = NCORES)
params2 <- attr(beta_lassosum2, "grid_param")

# --- 7. ALIGN WEIGHTS TO TARGET STUDY SAMPLE MATRIX (G) ---
message("Aligning lassosum2 weights with target study sample matrix (G)...")

# Clean df_beta of internal snp_match metadata before re-matching
df_beta_clean <- df_beta[, c("rsid", "chr", "pos", "a1", "a0", "beta", "beta_se", "n_eff")]

t_rsid <- tryCatch(snp_match(df_beta_clean, target_map, join_by_pos = FALSE), error = function(e) NULL)
t_pos  <- tryCatch(snp_match(df_beta_clean, target_map, join_by_pos = TRUE), error = function(e) NULL)

n_trsid <- if (!is.null(t_rsid)) nrow(t_rsid) else 0
n_tpos  <- if (!is.null(t_pos))  nrow(t_pos)  else 0

if (n_trsid >= n_tpos && n_trsid > 0) {
  target_match <- t_rsid
  message("Matched ", n_trsid, " variants with study sample matrix (G) by RSID.")
} else if (n_tpos > 0) {
  target_match <- t_pos
  message("Matched ", n_tpos, " variants with study sample matrix (G) by Position.")
} else {
  stop("Error: Failed to align lassosum2 variants with target study sample matrix G.")
}

idx_df <- target_match$`_NUM_ID_1`
idx_g  <- target_match$`_NUM_ID_`
if (is.null(idx_df)) {
  # bigsnpr >= 1.12 names the first-table index _NUM_ID_.ss (not _NUM_ID_1)
  idx_df <- target_match$`_NUM_ID_.ss`
}

if (is.null(idx_df) || is.null(idx_g) || length(idx_df) == 0) {
  stop("Error: Variant indexing failed during Step 7 matching.")
}

sign_flip <- ifelse(target_match$a1 == df_beta_clean$a1[idx_df], 1, -1)
sub_beta  <- beta_lassosum2[idx_df, , drop = FALSE]

final_beta   <- sweep(sub_beta, 1, sign_flip, "*")
final_g_cols <- idx_g

if (nrow(final_beta) != length(final_g_cols)) {
  stop(sprintf("Dimension mismatch before big_prodMat: final_beta has %d rows, final_g_cols has %d length.", 
               nrow(final_beta), length(final_g_cols)))
}

# --- 8. PREDICTION & TUNING ---
set.seed(args$seed)
n_total <- nrow(G)
valid_y_idx <- if (!is.null(y)) which(!is.na(y)) else 1:n_total

n_val <- min(args$n_val, floor(length(valid_y_idx) / 3))
if (n_val < 5) n_val <- min(length(valid_y_idx), 10)

ind.val  <- sample(valid_y_idx, n_val)
ind.test <- setdiff(rows_along(G), ind.val)

pred_grid2 <- big_prodMat(G, final_beta, ind.col = final_g_cols)

params2$score <- apply(pred_grid2[ind.val, , drop = FALSE], 2, function(x) {
  if (all(is.na(x)) || sd(x, na.rm = TRUE) == 0) return(NA_real_)
  if (is.null(y) || sd(y[ind.val], na.rm = TRUE) == 0) return(NA_real_)
  fit <- tryCatch(lm(y[ind.val] ~ x), error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  co <- summary(fit)$coefficients
  if (nrow(co) < 2) return(NA_real_)
  return(co[2, 3])
})

# --- 9. MODEL SELECTION & FALLBACK ---
best_idx <- which.max(params2$score)

if (length(best_idx) == 0 || is.na(best_idx) || all(is.na(params2$score))) {
  message("WARNING: Validation scores all NA. Choosing grid model with maximum prediction variance...")
  grid_vars <- apply(pred_grid2, 2, var, na.rm = TRUE)
  best_idx <- if (any(!is.na(grid_vars) & grid_vars > 0)) which.max(grid_vars) else ceiling(ncol(beta_lassosum2) / 2)
  best_col_name <- paste0("grid_", best_idx)
  test_r2 <- NA_real_
} else {
  best_col_name <- paste0("grid_", best_idx)
  pred_test <- pred_grid2[ind.test, best_idx]
  test_r2 <- tryCatch(pcor(pred_test, y[ind.test], NULL)^2, error = function(e) NA_real_)
}

best_beta <- final_beta[, best_idx]

p <- ggplot(params2, aes(x = lambda, y = score, color = as.factor(delta))) +
  theme_bigstatsr() + geom_point() + geom_line() + scale_x_log10() +
  labs(title = "lassosum2 Tuning", y = "Z-Score", color = "delta")
ggsave(paste0(args$out, "_lassosum_plot.png"), p)

write.csv(params2, paste0(args$out, "_grid_params.csv"), row.names = FALSE)
cat(paste("Best Test R2:", test_r2, "\n"), file = paste0(args$out, "_final_res.txt"))

unlink("temp_ld_lassosum", recursive = TRUE)

# --- 10. SAVE INDIVIDUAL SCORES & WEIGHTS ---
full_results <- data.frame(
  FID = obj.bigSNP$fam[, 1],
  IID = obj.bigSNP$fam[, 2],
  pred_grid2
)

colnames(full_results)[3:ncol(full_results)] <- paste0("grid_", 1:nrow(params2))
bigreadr::fwrite2(full_results, paste0(args$out, "_full_predictions.csv"))

final_prs_df <- full_results[, c("FID", "IID")]
final_prs_df$Best_PRS_Score <- full_results[[best_col_name]]
bigreadr::fwrite2(final_prs_df, paste0(args$out, "_final_best_prs.csv"))

export_map <- target_map[final_g_cols, ]

weights <- data.frame(
  SNP  = as.character(export_map$rsid),
  A1   = as.character(export_map$a1),
  BETA = best_beta,
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(SNP) & SNP != "" & SNP != "." &
         !is.na(A1) & A1 != "" &
         !is.na(BETA) & is.finite(BETA) & BETA != 0)

if (nrow(weights) == 0) {
  stop("Error: Weight table is empty after filtering non-finite/zero weights.")
}

weights_file <- paste0(args$out, "_weights.txt")
write.table(weights, weights_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

message(paste("Successfully exported weights to:", weights_file))

}, error = function(e) {
  message("\n========== ERROR IN run_lassosum2.R ==========")
  message("Condition: ", conditionMessage(e))
  message("Call stack:")
  for (i in seq_len(sys.nframe())) {
    call <- sys.call(i)
    if (!is.null(call)) message("  ", i, ": ", deparse(call)[1])
  }
  dump.frames("lassosum2_dump", to.file = TRUE)
  message("Dumped to lassosum2_dump.rda")
  stop(e)
})