#!/usr/bin/env Rscript

library(argparse)
library(bigsnpr)
library(ggplot2)
library(dplyr)
library(bigreadr)
library(data.table)

# 1. SET UP ARGUMENT PARSER
parser <- ArgumentParser(description='LDpred2 Pipeline for Polygenic Risk Scores')

parser$add_argument("--anc_bed", type="character", help="Path to the plink bed file")
parser$add_argument("--rds", type="character", help="Path to the .rds file")
parser$add_argument("--bim", type="character", required=TRUE, help="Path to the .bim file")

# Both are now optional; logic below handles various scenarios
parser$add_argument("--ss", type="character", required=TRUE, help="Path to summary statistics")
parser$add_argument("--beta_se", type="character", help="Optional: separate table with beta_se/SE")

parser$add_argument("--h2", type="numeric", default=0.4, help="Assumed heritability")
parser$add_argument("--n_val", type="integer", default=49, help="Number of samples for validation")
parser$add_argument("--out", type="character", default="ldpred2_out", help="Prefix for output files")

args <- parser$parse_args()

# --- 2. DATA LOADING & CONVERSION ---

if (!is.null(args$anc_bed)) {
  rds_path <- if (!is.null(args$rds)) args$rds else sub("\\.bed$", ".rds", args$anc_bed)
  if (!file.exists(rds_path)) {
    message("Converting .bed to .rds format...")
    snp_readBed(args$anc_bed)
  }
  args$rds <- rds_path 
}

if (is.null(args$rds)) stop("Error: You must provide either --anc_bed or --rds.")

obj.bigSNP <- snp_attach(args$rds)
G      <- obj.bigSNP$genotypes
CHR_id <- obj.bigSNP$map$chromosome
POS_id <- obj.bigSNP$map$physical.pos
y      <- obj.bigSNP$fam$affection
NCORES <- 1  # Avoid parallelism conflict with OpenMP

# Load BIM for alignment
bim.file <- fread2(args$bim, select = c(1, 4, 5, 6))
colnames(bim.file) <- c("chr", "pos", "bim.a1", "bim.a0")

# --- 3. SUMMARY STATISTICS PROCESSING ---

message("Reading summary statistics from: ", args$ss)
ss_raw <- fread2(args$ss)

# Map common GWAS headers to internal names
sumstats <- ss_raw %>%
  rename(
    rsid  = any_of(c("SNP", "rsid", "ID")),
    chr   = any_of(c("CHR", "chr", "chromosome")),
    a1    = any_of(c("A1", "a1", "allele1")),
    beta  = any_of(c("BETA", "beta", "eff")),
    n_eff = any_of(c("n_eff", "NMISS", "N"))
  )

# Parse position from SNP column (e.g., chr1:128161:G:C format)
if (!("pos" %in% colnames(sumstats)) && "rsid" %in% colnames(sumstats)) {
  message("Parsing position from SNP column...")
  sumstats$pos <- as.numeric(sapply(strsplit(sumstats$rsid, ":"), function(x) x[2]))
}

# Only filter by TEST if the column actually exists (e.g., PLINK outputs)
if ("TEST" %in% colnames(sumstats)) {
  sumstats <- sumstats %>% filter(TEST == "ADD")
}

# Logic for SE (Standard Error)
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

# Ensure required columns exist after mapping
req_cols <- c("chr", "pos", "a1", "beta", "beta_se", "n_eff")
missing <- setdiff(req_cols, colnames(sumstats))
if (length(missing) > 0) {
  stop("Missing required columns in summary stats: ", paste(missing, collapse = ", "))
}

# Merge with BIM to get a0 (reference allele)
sumstats <- sumstats %>%
  inner_join(bim.file, by = c("chr", "pos")) %>%
  mutate(a0 = ifelse(a1 == bim.a1, bim.a0, bim.a1)) %>%
  select(chr, pos, a1, a0, beta, beta_se, n_eff)

# --- 4. PREPARE FOR LDPRED2 ---

map <- setNames(obj.bigSNP$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
df_beta <- snp_match(sumstats, map, join_by_pos = TRUE)

POS2 <- obj.bigSNP$map$genetic.dist
ind.row <- rows_along(G)
maf <- snp_MAF(G, ind.row = ind.row, ind.col = df_beta$`_NUM_ID_`, ncores = NCORES)
maf_thr <- 1 / sqrt(length(ind.row))
df_beta <- df_beta[maf > maf_thr, ]

# --- 5. COMPUTE LD MATRIX ---

tmp <- tempfile(tmpdir = "temp_ld")
dir.create("temp_ld", showWarnings = FALSE)
corr <- NULL
ld <- NULL

for (chr in 1:22) {
  ind.chr <- which(df_beta$chr == chr)
  ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
  if (length(ind.chr2) < 2) next

  # Print progress for each chromosome
  message(paste("Processing Chromosome:", chr, "| SNPs:", length(ind.chr2)))
  
  corr0 <- snp_cor(G, ind.col = ind.chr2, size = 3/1000, infos.pos = POS2[ind.chr2], ncores = NCORES)

  if (is.null(corr)) {
    ld <- Matrix::colSums(corr0^2)
    corr <- as_SFBM(corr0, tmp, compact = TRUE)
  } else {
    ld <- c(ld, Matrix::colSums(corr0^2))
    corr$add_columns(corr0, nrow(corr))
  }
}

# --- 6. RUN MODELS ---

# Split validation/test sets
ind.val  <- sample(nrow(G), args$n_val)
ind.test <- setdiff(rows_along(G), ind.val)

# Infinitesimal Model
beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = args$h2)
pred_inf <- big_prodVec(G, beta_inf, ind.row = ind.test, ind.col = df_beta$`_NUM_ID_`)
r2_inf   <- pcor(pred_inf, y[ind.test], NULL)

# Grid Model
h2_seq <- round(args$h2 * c(0.7, 1, 1.4), 4)
p_seq  <- signif(seq_log(1e-4, 1, length.out = 10), 2)
params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE))

beta_grid <- snp_ldpred2_grid(corr, df_beta, params, ncores = NCORES)
pred_grid <- big_prodMat(G, beta_grid, ind.col = df_beta[["_NUM_ID_"]])

# Tuning on Validation set
params$score <- apply(pred_grid[ind.val, ], 2, function(x) {
  if (all(is.na(x))) return(NA)
  summary(lm(y[ind.val] ~ x))$coef[2, 3] # Use t-stat as score
})

# Save Grid Plot
p <- ggplot(params, aes(x = p, y = score, color = as.factor(h2))) +
  geom_point() + geom_line() + scale_x_log10() +
  facet_wrap(~ sparse) + labs(title="LDpred2 Grid Search Tuning")
ggsave(paste0(args$out, "_grid_plot.png"), p)

# Best Grid Prediction
best_grid_idx <- which.max(params$score)
pred_grid_best <- big_prodVec(G, beta_grid[, best_grid_idx], ind.row = ind.test, ind.col = df_beta[["_NUM_ID_"]])
r2_grid <- pcor(pred_grid_best, y[ind.test], NULL)

# --- 7. OUTPUT RESULTS ---

# 1. Calculate scores for EVERYONE (not just the test set)
# We already have pred_grid for everyone, we just need to pick the best column
prs_grid_all <- pred_grid[, best_grid_idx]

# Recalculate Infinitesimal for everyone (removing the ind.row restriction)
prs_inf_all <- big_prodVec(G, beta_inf, ind.col = df_beta[["_NUM_ID_"]])

# 2. Create the individual-level table
# obj.bigSNP$fam contains family.ID (FID) and sample.ID (IID)
prs_report <- data.frame(
  FID      = obj.bigSNP$fam$family.ID,
  IID      = obj.bigSNP$fam$sample.ID,
  PRS_inf  = prs_inf_all,
  PRS_grid = prs_grid_all
)

# 3. Save the individual scores
fwrite(prs_report, paste0(args$out, "_individual_scores.txt"), sep = "\t", row.names = FALSE)

# 4. Save the performance summary
results <- data.frame(
  Method = c("Infinitesimal", "Grid"),
  R2 = c(r2_inf^2, r2_grid^2)
)
write.csv(results, paste0(args$out, "_performance.csv"), row.names = FALSE)

message("Success! Individual scores saved to: ", paste0(args$out, "_individual_scores.txt"))