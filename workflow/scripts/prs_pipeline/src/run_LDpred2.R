#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(MKL_NUM_THREADS = "1")

library(argparse)
library(bigsnpr)
library(ggplot2)
library(dplyr)
library(bigreadr)
library(data.table)

options(bigstatsr.check.parallel = FALSE)
options(bigstatsr.check.args = FALSE)

# Helper for safe regression tuning on validation set
eval_score <- function(x, y_val) {
  if (all(is.na(x)) || sd(x, na.rm = TRUE) == 0) return(NA_real_)
  valid <- !is.na(y_val) & !is.na(x)
  if (sum(valid) < 5 || sd(y_val[valid]) == 0) return(NA_real_)
  fit <- tryCatch(lm(y_val[valid] ~ x[valid]), error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  co <- summary(fit)$coefficients
  if (nrow(co) < 2) return(NA_real_)
  return(co[2, 3]) # t-statistic
}

parser <- ArgumentParser(description='LDpred2 Pipeline for Polygenic Risk Scores')

parser$add_argument("--anc_bed", type="character", help="Path to the plink bed file")
parser$add_argument("--rds", type="character", help="Path to the .rds file")
parser$add_argument("--bim", type="character", required=TRUE, help="Path to the .bim file")
parser$add_argument("--ss", type="character", required=TRUE, help="Path to summary statistics")
parser$add_argument("--beta_se", type="character", help="Optional: separate table with beta_se/SE")
parser$add_argument("--afreq", type="character", help="Path to PLINK2 .afreq file")
parser$add_argument("--pheno", type="character", help="Path to phenotype file")

parser$add_argument("--h2", type="numeric", default=0.1, help="Fallback assumed heritability (default: 0.1)")
parser$add_argument("--n_val", type="integer", default=49, help="Number of samples for validation")
parser$add_argument("--seed", type="integer", default=1, help="Seed for validation/test split")
parser$add_argument("--out", type="character", default="ldpred2_out", help="Prefix for output files")
parser$add_argument("--ncores", type="integer", default=1, help="Number of CPU cores")
parser$add_argument("--ld-cache-dir", type="character", help="Directory to cache/reuse per-chromosome LD matrices")
parser$add_argument("--ld-matrix-dir", type="character", help="Directory with pre-computed LD matrix")

args <- parser$parse_args()

tryCatch({
  
  if (!is.null(args$anc_bed)) {
    rds_path <- if (!is.null(args$rds)) args$rds else sub("\\.bed$", ".rds", args$anc_bed)
    bk_path  <- sub("\\.bed$", ".bk", args$anc_bed)
    
    if (!file.exists(rds_path)) {
      if (file.exists(bk_path)) file.remove(bk_path)
      message("Converting .bed to .rds format...")
      snp_readBed(args$anc_bed)
    }
    args$rds <- rds_path 
  }
  
  if (is.null(args$rds)) stop("Error: You must provide either --anc_bed or --rds.")
  
  obj.bigSNP <- snp_attach(args$rds)
  G      <- obj.bigSNP$genotypes
  NCORES <- args$ncores
  y      <- obj.bigSNP$fam$affection
  
  if (!is.null(args$pheno)) {
    message("Reading phenotype from: ", args$pheno)
    pheno <- fread2(args$pheno)
    colnames(pheno)[1:3] <- c("FID", "IID", "pheno")
    
    fam <- obj.bigSNP$fam
    fam_key   <- paste(fam$family.ID, fam$sample.ID, sep = "_")
    pheno_key <- paste(pheno$FID, pheno$IID, sep = "_")
    
    m <- match(fam_key, pheno_key)
    if (any(is.na(m))) {
      warning("Missing phenotypes for ", sum(is.na(m)), " individuals in the genotype set.")
    }
    
    y <- pheno$pheno[m]
    message("Successfully matched ", sum(!is.na(y)), " / ", length(y), " individuals to phenotype data.")
  }
  
  message("Reading summary statistics from: ", args$ss)
  ss_raw <- fread2(args$ss)
  
  sumstats <- ss_raw %>%
    rename(
      rsid  = any_of(c("SNP", "rsid", "ID")),
      pos   = any_of(c("BP", "pos", "position")),
      chr   = any_of(c("CHR", "chr", "chromosome")),
      a1    = any_of(c("A1", "a1", "allele1")),
      a0    = any_of(c("A2", "a0", "allele2", "REF")),
      beta  = any_of(c("BETA", "beta", "eff")),
      n_eff = any_of(c("n_eff", "NMISS", "N"))
    )
  
  if ("TEST" %in% colnames(sumstats)) {
    sumstats <- sumstats %>% filter(TEST == "ADD")
  }
  
  if (!is.null(args$beta_se)) {
    message("Merging with external SE file: ", args$beta_se)
    beta_se_tab <- fread2(args$beta_se) %>% 
      rename(rsid = any_of(c("SNP", "rsid")), 
             beta_se = any_of(c("SE", "beta_se", "beta_se_tab"))) %>%
      select(rsid, beta_se)
    
    sumstats <- sumstats %>% inner_join(beta_se_tab, by = "rsid")
  } else {
    message("Extracting SE from primary summary stat file...")
    sumstats <- sumstats %>% rename(beta_se = any_of(c("SE", "beta_se")))
  }
  
  req_cols <- c("chr", "pos", "a1", "a0", "beta", "beta_se", "n_eff")
  missing <- setdiff(req_cols, colnames(sumstats))
  if (length(missing) > 0) {
    stop("Missing required columns in summary stats: ", paste(missing, collapse = ", "))
  }
  
  if (!is.null(args$ld_matrix_dir)) {
    message("Aligning sumstats with pre-computed LD matrix map...")
    map_ldref <- readRDS(file.path(args$ld_matrix_dir, "map.rds"))
  } else {
    map_ldref <- setNames(obj.bigSNP$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
  }
  
  # Step 1: Alignment via snp_match against LD reference
  info_snp <- snp_match(sumstats, map_ldref, join_by_pos = FALSE)
  info_snp <- tidyr::drop_na(tibble::as_tibble(info_snp))
  
  # Step 2: Quality Control - Auto-scaled SD Filter
  if (!is.null(args$ld_matrix_dir) && "af_UKBB" %in% names(info_snp)) {
    message("Applying SD QC filter between summary stats and LD reference...")
    sd_ldref  <- with(info_snp, sqrt(2 * af_UKBB * (1 - af_UKBB)))
    sd_ss_raw <- with(info_snp, 2 / sqrt(n_eff * beta_se^2))
    
    scale_factor <- median(sd_ldref / sd_ss_raw, na.rm = TRUE)
    if (!is.finite(scale_factor) || scale_factor <= 0) scale_factor <- 1
    sd_ss <- sd_ss_raw * scale_factor
    
    is_bad <- sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.05 | sd_ldref < 0.01
    is_bad[is.na(is_bad)] <- TRUE
    
    bad_count <- sum(is_bad)
    if (bad_count > (0.5 * nrow(info_snp))) {
      message(sprintf("WARNING: SD QC filter flagged %d / %d variants (>50%%). Skipping SD QC filter to prevent variant loss.", bad_count, nrow(info_snp)))
      df_beta <- info_snp
    } else {
      message(sprintf("Filtering out %d variants due to SD misalignment QC.", bad_count))
      df_beta <- info_snp[!is_bad, ]
    }
  } else {
    df_beta <- info_snp
  }
  
  # Step 3: Restrict to variants present in target dataset G
  target_map <- setNames(obj.bigSNP$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
  in_target  <- vctrs::vec_in(df_beta[, c("chr", "pos")], target_map[, c("chr", "pos")])
  df_beta    <- df_beta[in_target, ]
  
  CHRS <- 1:22
  df_beta <- df_beta[df_beta$chr %in% CHRS, ]
  message("Autosomes 1-22 after filtering: ", nrow(df_beta))
  
  if (nrow(df_beta) < 10) {
    stop("Fewer than 10 variants remain after filtering.")
  }
  
  # Step 4: Construct LD Matrix (SFBM)
  tmp <- tempfile(tmpdir = "temp_ld")
  dir.create("temp_ld", showWarnings = FALSE)
  corr <- NULL
  keep_idx <- logical(nrow(df_beta))
  
  if (!is.null(args$ld_matrix_dir)) {
    message("Loading pre-computed LD matrix from: ", args$ld_matrix_dir)
    ld_map <- readRDS(file.path(args$ld_matrix_dir, "map.rds"))
    
    for (chr in CHRS) {
      ind.chr <- which(df_beta$chr == chr)
      if (length(ind.chr) < 2) next
      
      map_chr_idx <- which(ld_map$chr == chr)
      local_idx   <- match(df_beta$`_NUM_ID_`[ind.chr], map_chr_idx)
      
      bad <- which(is.na(local_idx))
      if (length(bad) > 0) {
        local_idx <- local_idx[-bad]
        ind.chr   <- ind.chr[-bad]
      }
      if (length(local_idx) < 2) next
      
      keep_idx[ind.chr] <- TRUE
      corr0_full <- readRDS(file.path(args$ld_matrix_dir, sprintf("chr%d_corr.rds", chr)))
      corr0      <- corr0_full[local_idx, local_idx]
      
      if (any(is.na(corr0))) corr0[is.na(corr0)] <- 0
      corr0 <- corr0 + diag(ncol(corr0)) * 1e-5
      
      if (is.null(corr)) {
        corr <- as_SFBM(corr0, tmp, compact = TRUE)
      } else {
        corr$add_columns(corr0, nrow(corr))
      }
    }
  } else {
    POS2 <- if ("genetic.dist" %in% names(obj.bigSNP$map) && !all(is.na(obj.bigSNP$map$genetic.dist))) {
      obj.bigSNP$map$genetic.dist
    } else {
      obj.bigSNP$map$physical.pos / 1e6
    }
    
    for (chr in CHRS) {
      ind.chr  <- which(df_beta$chr == chr)
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
    }
  }
  
  df_beta <- df_beta[keep_idx, ]
  message("Variants in LD matrix: ", nrow(df_beta))
  
  if (is.null(corr) || nrow(df_beta) != ncol(corr)) {
    stop("Mismatch between df_beta rows and corr cols after LD construction.")
  }
  
  # Step 5: Estimate Heritability via LDSC
  h2_est <- args$h2
  if ("ld" %in% names(df_beta)) {
    message("Estimating heritability with snp_ldsc...")
    ldsc <- tryCatch({
      snp_ldsc(df_beta$ld, 
               ld_size = nrow(map_ldref),
               chi2 = (df_beta$beta / df_beta$beta_se)^2,
               sample_size = df_beta$n_eff, 
               ncores = NCORES)
    }, error = function(e) NULL)
    
    if (!is.null(ldsc) && !is.na(ldsc[["h2"]]) && ldsc[["h2"]] > 0) {
      h2_est <- ldsc[["h2"]]
      message(sprintf("LDSC estimated h2: %.4f", h2_est))
    } else {
      message(sprintf("LDSC estimation failed or returned non-positive h2. Using fallback h2: %.4f", h2_est))
    }
  }

  # --- Run LDpred2 Models ---
  
  # 1. Infinitesimal Model
  message("Running LDpred2-inf...")
  beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
  
  # 2. Grid Model
  message("Running LDpred2-grid...")
  h2_seq <- pmax(0.0001, pmin(0.95, round(h2_est * c(0.1, 0.3, 0.7, 1, 1.4), 4)))
  p_seq  <- signif(seq_log(1e-4, 1, length.out = 10), 2)
  params <- expand.grid(p = p_seq, h2 = unique(h2_seq), sparse = c(FALSE, TRUE))
  
  beta_grid <- snp_ldpred2_grid(corr, df_beta, params, ncores = NCORES)

  # Step 6: Final Alignment to Target Dataset
  # Strip ALL prior internal bigsnpr columns and attach explicit row counter
  df_beta_clean <- df_beta %>% 
    transmute(
      chr       = chr,
      pos       = pos,
      a0        = a0,
      a1        = a1,
      beta      = beta,
      beta_se   = beta_se,
      n_eff     = n_eff,
      rsid      = rsid,
      ldpred_id = seq_len(n())
    )
  
  df_target <- snp_match(df_beta_clean, target_map, join_by_pos = FALSE)
  target_id <- df_target$`_NUM_ID_`      # Column index in target matrix G
  local_id  <- df_target$ldpred_id       # Row index in df_beta / beta_grid (1..nrow(df_beta))
  flip      <- df_target$beta            # Allele flip indicator (+1 or -1)
  
  beta_inf_target  <- beta_inf[local_id] * flip
  beta_grid_target <- beta_grid[local_id, , drop = FALSE] * flip
  
  set.seed(args$seed)
  n_total <- nrow(G)
  valid_pheno_idx <- which(!is.na(y))
  if (length(valid_pheno_idx) < 10) stop("Too few samples with non-missing phenotype data.")
  
  n_val    <- min(args$n_val, floor(length(valid_pheno_idx) / 3))
  ind.val  <- sample(valid_pheno_idx, n_val)
  ind.test <- setdiff(rows_along(G), ind.val)
  
  # Evaluation & Predictions
  pred_inf  <- big_prodVec(G, beta_inf_target, ind.row = ind.test, ind.col = target_id)
  r2_inf    <- pcor(pred_inf, y[ind.test], NULL)
  
  pred_grid <- big_prodMat(G, beta_grid_target, ind.col = target_id)
  
  # Tune Grid parameters on validation set
  params$score <- apply(pred_grid[ind.val, , drop = FALSE], 2, function(x) eval_score(x, y[ind.val]))
  
  p <- ggplot(params, aes(x = p, y = score, color = as.factor(h2))) +
    geom_point() + geom_line() + scale_x_log10() +
    facet_wrap(~ sparse) + labs(title = "LDpred2 Grid Search Tuning")
  ggsave(paste0(args$out, "_grid_plot.png"), p)
  
  best_grid_idx <- which.max(params$score)
  if (length(best_grid_idx) == 0 || is.na(best_grid_idx) || all(is.na(params$score))) {
    message("WARNING: All grid model scores are NA. Using default/fallback grid index.")
    best_grid_idx  <- 1
    pred_grid_best <- rep(NA_real_, length(ind.test))
    r2_grid        <- NA_real_
    prs_grid_all   <- rep(NA_real_, n_total)
  } else {
    pred_grid_best <- pred_grid[ind.test, best_grid_idx]
    r2_grid        <- pcor(pred_grid_best, y[ind.test], NULL)
    prs_grid_all   <- pred_grid[, best_grid_idx]
  }
  
  # --- Output Results ---
  prs_inf_all <- big_prodVec(G, beta_inf_target, ind.col = target_id)
  
  prs_report <- data.frame(
    FID      = obj.bigSNP$fam$family.ID,
    IID      = obj.bigSNP$fam$sample.ID,
    PRS_inf  = prs_inf_all,
    PRS_grid = prs_grid_all
  )
  fwrite(prs_report, paste0(args$out, "_individual_scores.txt"), sep = "\t", row.names = FALSE)
  
  results <- data.frame(
    Method = c("Infinitesimal", "Grid"),
    R2 = c(r2_inf^2, r2_grid^2)
  )
  write.csv(results, paste0(args$out, "_performance.csv"), row.names = FALSE)
  
  map_sub <- obj.bigSNP$map[target_id, ]
  
  inf_weights <- data.frame(
    SNP  = as.character(map_sub$marker.ID),
    A1   = as.character(map_sub$allele1),
    BETA = beta_inf_target,
    stringsAsFactors = FALSE
  )
  fwrite(inf_weights, paste0(args$out, "_inf_weights.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  
  grid_weights <- data.frame(
    SNP  = as.character(map_sub$marker.ID),
    A1   = as.character(map_sub$allele1),
    BETA = beta_grid_target[, best_grid_idx],
    stringsAsFactors = FALSE
  )
  fwrite(grid_weights, paste0(args$out, "_grid_weights.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  
  message("Success! Individual scores saved to: ", paste0(args$out, "_individual_scores.txt"))  
  
}, error = function(e) {
  message("\n========== ERROR IN run_LDpred2.R ==========")
  message("Condition: ", conditionMessage(e))
  dump.frames("ldpred2_dump", to.file = TRUE)
  stop(e)
})