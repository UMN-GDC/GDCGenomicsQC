#!/usr/bin/env Rscript

library(argparse)
library(tidyverse)
library(mvtnorm)

parser <- ArgumentParser(description = "Run complete phenotype simulation pipeline")
parser$add_argument("--ancestry1", type = "character", required = TRUE,
    help = "Filepath to plink files for ancestry 1 (base, no extension)")
parser$add_argument("--ancestry2", type = "character", required = TRUE,
    help = "Filepath to plink files for ancestry 2 (base, no extension)")
parser$add_argument("--out_dir", type = "character", required = TRUE,
    help = "Output directory for simulation files")
parser$add_argument("--anc1_name", type = "character", required = TRUE,
    help = "Ancestry 1 name (e.g., AFR, EUR)")
parser$add_argument("--anc2_name", type = "character", required = TRUE,
    help = "Ancestry 2 name (e.g., AFR, EUR)")
parser$add_argument("--n_sims", type = "integer", default = 10,
    help = "Number of simulations to run")
parser$add_argument("--seed", type = "integer", default = 42,
    help = "Random seed")
parser$add_argument("--heritability", type = "numeric", default = 0.4,
    help = "Heritability (default 0.4)")
parser$add_argument("--rho", type = "numeric", default = 0.8,
    help = "Genetic correlation (default 0.8)")
parser$add_argument("--maf", type = "numeric", default = 0.05,
    help = "Minor allele frequency cutoff (default 0.05)")
parser$add_argument("--skip_thinning", type = "character", default = "true",
    help = "Skip thinning (default true)")
parser$add_argument("--thin_count_snps", type = "integer", default = 1000000,
    help = "Number of SNPs to thin to")
parser$add_argument("--thin_count_inds", type = "integer", default = 10000,
    help = "Number of individuals to thin to")

args <- parser$parse_args()

cat("========================================\n")
cat("Phenotype Simulation Pipeline\n")
cat("========================================\n")
cat("Ancestry 1:", args$anc1_name, "\n")
cat("Ancestry 2:", args$anc2_name, "\n")
cat("Number of simulations:", args$n_sims, "\n")
cat("Heritability:", args$heritability, "\n")
cat("Genetic correlation:", args$rho, "\n")
cat("MAF threshold:", args$maf, "\n")
cat("Seed:", args$seed, "\n")
cat("Skip thinning:", args$skip_thinning, "\n")
cat("Output directory:", args$out_dir, "\n")
cat("========================================\n\n")

dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)

anc1File <- args$ancestry1
anc2File <- args$ancestry2
herit <- args$heritability
rho <- args$rho
maf <- args$maf
n_sims <- args$n_sims
seed <- args$seed

cat("STEP 1: Thinning and preprocessing\n")
cat("-----------------------------------\n")

skip_thinning <- tolower(args$skip_thinning) == "true" | tolower(args$skip_thinning) == "1"

if (skip_thinning) {
    cat("Skipping thinning step\n")
    system2("plink2", args = c("--pfile", anc1File, "--maf", as.character(maf), "--freq", "--out", anc1File))
    system2("awk", args = c("{print $2}", paste0(anc1File, ".pvar")), stdout = paste0(anc1File, ".snpIDs"))
    system2("plink2", args = c("--pfile", anc2File, "--maf", as.character(maf), "--freq", "--out", anc2File))
} else {
    cat("Applying thinning\n")
    system2("plink2", args = c("--pfile", anc1File, "--thin-count", as.character(args$thin_count_snps),
        "--thin-indiv-count", as.character(args$thin_count_inds), "--maf", as.character(maf),
        "--freq", "--out", anc1File))
    system2("awk", args = c("{print $2}", paste0(anc1File, ".pvar")), stdout = paste0(anc1File, ".snpIDs"))
    system2("plink2", args = c("--pfile", anc2File, "--extract", paste0(anc1File, ".snpIDs"),
        "--thin-indiv-count", as.character(args$thin_count_inds), "--maf", as.character(maf),
        "--freq", "--out", anc2File))
}

cat("\nSTEP 2: Creating base .fam files\n")
cat("-----------------------------------\n")

system2("plink2", args = c("--pfile", anc1File, "--keep-allele-order", "--make-bed",
    "--out", paste0(anc1File, "_base")))
system2("plink2", args = c("--pfile", anc2File, "--keep-allele-order", "--make-bed",
    "--out", paste0(anc2File, "_base")))

base_fam_1 <- paste0(args$out_dir, "/", args$anc1_name, "_base.fam")
base_fam_2 <- paste0(args$out_dir, "/", args$anc2_name, "_base.fam")

if (!file.exists(paste0(anc1File, "_base.fam"))) {
    stop(paste0("Base .fam file not found: ", anc1File, "_base.fam"))
}
if (!file.exists(paste0(anc2File, "_base.fam"))) {
    stop(paste0("Base .fam file not found: ", anc2File, "_base.fam"))
}

# Create base .fam files with only the first 5 columns (FID, IID, PID, MID, Sex)
cat("Trimming .fam files to first 5 columns...\n")
fam1_base <- read_table(paste0(anc1File, "_base.fam"), col_names = FALSE)
fam2_base <- read_table(paste0(anc2File, "_base.fam"), col_names = FALSE)

write_delim(fam1_base[, 1:5], base_fam_1, col_names = FALSE, delim = " ")
write_delim(fam2_base[, 1:5], base_fam_2, col_names = FALSE, delim = " ")

sim_fam_1 <- paste0(args$out_dir, "/", args$anc1_name, "_simulation.fam")
sim_fam_2 <- paste0(args$out_dir, "/", args$anc2_name, "_simulation.fam")

file.copy(base_fam_1, sim_fam_1, overwrite = TRUE)
file.copy(base_fam_2, sim_fam_2, overwrite = TRUE)

cat("Base .fam files created\n\n")

simEffects <- function(anc1_afreq, anc2_afreq, herit, rho, maf) {
    anc1 <- read_table(anc1_afreq)
    C_ANC1 <- sum(anc1$ALT_FREQS > maf, na.rm = TRUE)
    anc2 <- read_table(anc2_afreq)
    C_ANC2 <- sum(anc2$ALT_FREQS > maf, na.rm = TRUE)

    combined <- inner_join(anc1, anc2, by = c("#CHROM", "ID", "REF", "ALT"), suffix = c("_1", "_2"))

    combined <- combined %>%
        mutate(
            sdEff_1 = case_when(
                is.na(ALT_FREQS_1) ~ 0,
                ALT_FREQS_1 < maf ~ 0,
                .default = sqrt(herit / C_ANC1)
            ),
            sdEff_2 = case_when(
                is.na(ALT_FREQS_2) ~ 0,
                ALT_FREQS_2 < maf ~ 0,
                .default = sqrt(herit / C_ANC2)
            )
        )

    betas <- matrix(NA, nrow = nrow(combined), ncol = 2)
    cov_matrix <- matrix(c(
        combined$sdEff_1[1], (rho * herit) / (sqrt(C_ANC1) * sqrt(C_ANC2)),
        (rho * herit) / (sqrt(C_ANC2) * sqrt(C_ANC1)), combined$sdEff_2[1]
    ), ncol = 2)

    for (i in 1:nrow(combined)) {
        cov_matrix[1, 1] <- combined$sdEff_1[i]
        cov_matrix[2, 2] <- combined$sdEff_2[i]
        cov_matrix[1, 2] <- (rho * herit) / (sqrt(C_ANC1) * sqrt(C_ANC2))
        cov_matrix[2, 1] <- (rho * herit) / (sqrt(C_ANC2) * sqrt(C_ANC1))

        if (any(is.na(cov_matrix)) || any(is.nan(cov_matrix))) {
            betas[i, ] <- c(0, 0)
        } else {
            betas[i, 1:2] <- rmvnorm(n = 1, mean = c(0, 0), sigma = cov_matrix)
        }
    }

    combined$effs_1 <- betas[, 1] * sqrt(2 * combined$ALT_FREQS_1 * (1 - combined$ALT_FREQS_1))
    combined$effs_2 <- betas[, 2] * sqrt(2 * combined$ALT_FREQS_2 * (1 - combined$ALT_FREQS_2))

    combined$effs_1 <- combined$effs_1 * sqrt(herit / sum(combined$effs_1^2, na.rm = TRUE))
    combined$effs_2 <- combined$effs_2 * sqrt(herit / sum(combined$effs_2^2, na.rm = TRUE))

    return(combined)
}

adjustPheno <- function(effects_file, sscore_file, herit) {
    pheno <- read_table(sscore_file)
    sd_score <- sd(pheno$SCORE1_AVG, na.rm = TRUE)

    effects <- read_delim(effects_file, delim = "\t")

    if (sd_score > 0) {
        effects[, 3] <- effects[, 3] / sd_score * sqrt(herit)
    }

    write_delim(effects, effects_file, delim = "\t")
}

finalizePheno <- function(sscore_file, pheno_file, herit) {
    pheno <- read_table(sscore_file)
    pheno <- pheno %>%
        mutate(
            error = rnorm(n = nrow(pheno), 0, sqrt(1 - herit)),
            error = error / sd(error, na.rm = TRUE) * sqrt(1 - herit),
            pheno = SCORE1_AVG + error
        ) %>%
        select(-error)

    write_delim(pheno, pheno_file, delim = "\t")
}

cat("STEP 3: Running simulations\n")
cat("-----------------------------------\n")

for (sim_i in 1:n_sims) {
    cat(sprintf("Running simulation %d of %d...\n", sim_i, n_sims))

    current_seed <- seed + sim_i - 1
    set.seed(current_seed)

    cat(sprintf("  [Sim %d] Generating SNP effect sizes...\n", sim_i))
    combined_effs <- simEffects(paste0(anc1File, ".afreq"), paste0(anc2File, ".afreq"), herit, rho, maf)

    write_delim(combined_effs %>% select(`#ID` = ID, REF, effs_1), paste0(anc1File, "Effects.csv"), delim = "\t")
    write_delim(combined_effs %>% select(`#ID` = ID, REF, effs_2), paste0(anc2File, "Effects.csv"), delim = "\t")

    cat(sprintf("  [Sim %d] Computing PRS (round 1)...\n", sim_i))
    system2("plink2", args = c("--pfile", anc1File, "--score", paste0(anc1File, "Effects.csv"),
        "1", "2", "3", "header", "--out", anc1File))
    system2("plink2", args = c("--pfile", anc2File, "--score", paste0(anc2File, "Effects.csv"),
        "1", "2", "3", "header", "--out", anc2File))

    cat(sprintf("  [Sim %d] Adjusting phenotype for heritability...\n", sim_i))
    adjustPheno(paste0(anc1File, "Effects.csv"), paste0(anc1File, ".sscore"), herit)
    adjustPheno(paste0(anc2File, "Effects.csv"), paste0(anc2File, ".sscore"), herit)

    cat(sprintf("  [Sim %d] Computing PRS (round 2)...\n", sim_i))
    system2("plink2", args = c("--pfile", anc1File, "--score", paste0(anc1File, "Effects.csv"),
        "1", "2", "3", "header", "--out", anc1File))
    system2("plink2", args = c("--pfile", anc2File, "--score", paste0(anc2File, "Effects.csv"),
        "1", "2", "3", "header", "--out", anc2File))

    cat(sprintf("  [Sim %d] Finalizing phenotype with noise...\n", sim_i))
    finalizePheno(paste0(anc1File, ".sscore"), paste0(anc1File, ".pheno"), herit)
    finalizePheno(paste0(anc2File, ".sscore"), paste0(anc2File, ".pheno"), herit)

    # Read the last column from the .pheno file
    pheno_1_data <- read_table(paste0(anc1File, ".pheno"), col_names = TRUE)
    pheno_2_data <- read_table(paste0(anc2File, ".pheno"), col_names = TRUE)
    
    pheno_1 <- pheno_1_data[[ncol(pheno_1_data)]]
    pheno_2 <- pheno_2_data[[ncol(pheno_2_data)]]

    # Read the current .fam file
    fam_1 <- read_table(sim_fam_1, col_names = FALSE)
    fam_2 <- read_table(sim_fam_2, col_names = FALSE)

    # Add the new phenotype column
    fam_1 <- fam_1 %>% mutate(!!paste0("pheno_", sim_i) := pheno_1)
    fam_2 <- fam_2 %>% mutate(!!paste0("pheno_", sim_i) := pheno_2)

    # Write the updated .fam file (space-separated, as is typical for PLINK)
    write_delim(fam_1, sim_fam_1, col_names = FALSE, delim = " ")
    write_delim(fam_2, sim_fam_2, col_names = FALSE, delim = " ")

    cat(sprintf("  [Sim %d] Phenotype appended to .fam file (Total columns: %d)\n", sim_i, ncol(fam_1)))
}

cat("\nSTEP 4: Creating PLINK binary files\n")
cat("-----------------------------------\n")

# Save the multi-column .fam files to temp names
temp_fam_1 <- paste0(sim_fam_1, ".tmp")
temp_fam_2 <- paste0(sim_fam_2, ".tmp")
file.copy(sim_fam_1, temp_fam_1, overwrite = TRUE)
file.copy(sim_fam_2, temp_fam_2, overwrite = TRUE)

cat("Creating .bed/.bim for ancestry 1...\n")
system2("plink2", args = c("--pfile", anc1File, "--keep-allele-order", "--make-bed",
    "--out", paste0(args$out_dir, "/", args$anc1_name, "_simulation")))

cat("Creating .bed/.bim for ancestry 2...\n")
system2("plink2", args = c("--pfile", anc2File, "--keep-allele-order", "--make-bed",
    "--out", paste0(args$out_dir, "/", args$anc2_name, "_simulation")))

# Restore the multi-column .fam files
file.rename(temp_fam_1, sim_fam_1)
file.rename(temp_fam_2, sim_fam_2)

cat("\n========================================\n")
cat("Simulation complete!\n")
cat("========================================\n")

fam_1_final <- read_table(sim_fam_1, col_names = FALSE)
fam_2_final <- read_table(sim_fam_2, col_names = FALSE)
n_pheno_cols_1 <- ncol(fam_1_final) - 5
n_pheno_cols_2 <- ncol(fam_2_final) - 5

cat(sprintf("Output files:\n"))
cat(sprintf("  %s: %d samples, %d phenotype columns\n", sim_fam_1, nrow(fam_1_final), n_pheno_cols_1))
cat(sprintf("  %s: %d samples, %d phenotype columns\n", sim_fam_2, nrow(fam_2_final), n_pheno_cols_2))