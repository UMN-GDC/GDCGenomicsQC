#!/usr/bin/env Rscript

library(argparse)
library(tidyverse)

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

validate_simulation_parameters <- function(herit, rho, maf) {
    if (herit <= 0 || herit >= 1) stop("heritability must be strictly between 0 and 1")
    if (rho < -1 || rho > 1) stop("rho must be between -1 and 1")
    if (maf <= 0 || maf >= 0.5) stop("maf must be strictly between 0 and 0.5")
}

run_or_stop <- function(command, command_args) {
    status <- system2(command, args = command_args)
    if (!identical(status, 0L)) {
        stop(sprintf("Command failed (%s): %s", status, paste(c(command, command_args), collapse = " ")))
    }
}

get_score_ids <- function(score_data) {
    fid_col <- intersect(c("#FID", "FID", "X.FID"), names(score_data))
    if (length(fid_col) != 1 || !all(c("IID", "SCORE1_AVG") %in% names(score_data))) {
        stop("PLINK score file must contain FID, IID, and SCORE1_AVG columns")
    }
    fid_col <- fid_col[[1]]
    tibble(
        FID = as.character(score_data[[fid_col]]),
        IID = as.character(score_data[["IID"]]),
        genetic_score = score_data[["SCORE1_AVG"]]
    )
}

simEffects <- function(anc1_afreq, anc2_afreq, herit, rho, maf) {
    anc1 <- read_table(anc1_afreq, show_col_types = FALSE)
    anc2 <- read_table(anc2_afreq, show_col_types = FALSE)
    combined <- inner_join(anc1, anc2, by = c("#CHROM", "ID", "REF", "ALT"), suffix = c("_1", "_2"))

    if (nrow(combined) == 0) stop("No variants overlap between ancestry frequency files")

    combined <- combined %>%
        mutate(
            maf_1 = pmin(ALT_FREQS_1, 1 - ALT_FREQS_1),
            maf_2 = pmin(ALT_FREQS_2, 1 - ALT_FREQS_2),
            eligible_1 = !is.na(maf_1) & maf_1 >= maf,
            eligible_2 = !is.na(maf_2) & maf_2 >= maf,
            shared_eligible = eligible_1 & eligible_2
        )

    c_anc1 <- sum(combined$eligible_1)
    c_anc2 <- sum(combined$eligible_2)
    if (c_anc1 == 0 || c_anc2 == 0) stop("No variants pass the MAF threshold in one or both ancestries")
    if (sum(combined$shared_eligible) < 2) stop("At least two shared MAF-eligible variants are required to simulate rho")

    # Effects are generated on the standardized-genotype scale. Shared eligible
    # variants use correlated normal draws; ancestry-specific variants are independent.
    z1 <- rnorm(nrow(combined))
    z2 <- rnorm(nrow(combined))
    beta_std_1 <- ifelse(combined$eligible_1, sqrt(herit / c_anc1) * z1, 0)
    beta_std_2 <- ifelse(combined$eligible_2, sqrt(herit / c_anc2) * z2, 0)
    beta_std_2[combined$shared_eligible] <- sqrt(herit / c_anc2) * (
        rho * z1[combined$shared_eligible] + sqrt(1 - rho^2) * z2[combined$shared_eligible]
    )

    combined$beta_std_1 <- beta_std_1
    combined$beta_std_2 <- beta_std_2
    combined$effs_1 <- ifelse(combined$eligible_1, beta_std_1 / sqrt(2 * combined$ALT_FREQS_1 * (1 - combined$ALT_FREQS_1)), 0)
    combined$effs_2 <- ifelse(combined$eligible_2, beta_std_2 / sqrt(2 * combined$ALT_FREQS_2 * (1 - combined$ALT_FREQS_2)), 0)
    combined
}

write_score_weights <- function(effects, ancestry, effects_file) {
    effect_column <- if (ancestry == 1) "effs_1" else "effs_2"
    weights <- effects %>%
        transmute(`#ID` = ID, REF, BETA = .data[[effect_column]])
    write_delim(weights, effects_file, delim = "\t")
}

calibrate_score_weights <- function(effects_file, sscore_file, herit) {
    score_data <- read_table(sscore_file, show_col_types = FALSE)
    score_sd <- sd(score_data$SCORE1_AVG, na.rm = TRUE)
    if (!is.finite(score_sd) || score_sd <= 0) stop("Initial genetic score has zero or non-finite variance")

    effects <- read_delim(effects_file, delim = "\t", show_col_types = FALSE)
    effects$BETA <- effects$BETA / score_sd * sqrt(herit)
    write_delim(effects, effects_file, delim = "\t")
}

finalizePheno <- function(sscore_file, herit) {
    genetic_score <- get_score_ids(read_table(sscore_file, show_col_types = FALSE))
    residual <- rnorm(nrow(genetic_score), 0, sqrt(1 - herit))
    residual <- residual / sd(residual) * sqrt(1 - herit)

    genetic_score %>%
        mutate(
            residual = residual,
            phenotype = genetic_score + residual,
            realized_h2 = var(genetic_score) / var(phenotype)
        )
}

append_phenotype_to_fam <- function(sim_fam, phenotype_data, sim_i) {
    fam <- read_table(sim_fam, col_names = FALSE, show_col_types = FALSE)
    fam_ids <- transmute(fam, FID = as.character(X1), IID = as.character(X2))
    matched <- match(paste(fam_ids$FID, fam_ids$IID), paste(phenotype_data$FID, phenotype_data$IID))
    if (anyNA(matched)) stop(sprintf("Could not match every individual in %s to its simulated phenotype", sim_fam))

    fam <- fam %>% mutate(!!paste0("pheno_", sim_i) := phenotype_data$phenotype[matched])
    write_delim(fam, sim_fam, col_names = FALSE, delim = " ")
    fam
}

cat("STEP 3: Running simulations\n")
cat("-----------------------------------\n")

validate_simulation_parameters(herit, rho, maf)
truth_dir <- file.path(args$out_dir, "simulation_truth")
work_dir <- file.path(args$out_dir, "simulation_work")
dir.create(truth_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
simulation_manifest <- vector("list", n_sims)

for (sim_i in 1:n_sims) {
    cat(sprintf("Running simulation %d of %d...\n", sim_i, n_sims))

    current_seed <- seed + sim_i - 1
    set.seed(current_seed)

    cat(sprintf("  [Sim %d] Generating SNP effect sizes...\n", sim_i))
    combined_effs <- simEffects(paste0(anc1File, ".afreq"), paste0(anc2File, ".afreq"), herit, rho, maf)

    effects_file_1 <- file.path(work_dir, paste0(args$anc1_name, "_score_weights.tsv"))
    effects_file_2 <- file.path(work_dir, paste0(args$anc2_name, "_score_weights.tsv"))
    score_prefix_1 <- file.path(work_dir, paste0(args$anc1_name, "_sim", sprintf("%03d", sim_i), "_score"))
    score_prefix_2 <- file.path(work_dir, paste0(args$anc2_name, "_sim", sprintf("%03d", sim_i), "_score"))
    write_score_weights(combined_effs, 1, effects_file_1)
    write_score_weights(combined_effs, 2, effects_file_2)

    cat(sprintf("  [Sim %d] Computing PRS (round 1)...\n", sim_i))
    run_or_stop("plink2", c("--pfile", anc1File, "--score", effects_file_1, "1", "2", "3", "header", "--out", score_prefix_1))
    run_or_stop("plink2", c("--pfile", anc2File, "--score", effects_file_2, "1", "2", "3", "header", "--out", score_prefix_2))

    cat(sprintf("  [Sim %d] Adjusting phenotype for heritability...\n", sim_i))
    calibrate_score_weights(effects_file_1, paste0(score_prefix_1, ".sscore"), herit)
    calibrate_score_weights(effects_file_2, paste0(score_prefix_2, ".sscore"), herit)

    cat(sprintf("  [Sim %d] Computing PRS (round 2)...\n", sim_i))
    run_or_stop("plink2", c("--pfile", anc1File, "--score", effects_file_1, "1", "2", "3", "header", "--out", score_prefix_1))
    run_or_stop("plink2", c("--pfile", anc2File, "--score", effects_file_2, "1", "2", "3", "header", "--out", score_prefix_2))

    cat(sprintf("  [Sim %d] Finalizing phenotype with noise...\n", sim_i))
    pheno_1_data <- finalizePheno(paste0(score_prefix_1, ".sscore"), herit)
    pheno_2_data <- finalizePheno(paste0(score_prefix_2, ".sscore"), herit)
    write_delim(pheno_1_data, file.path(truth_dir, sprintf("%s_sim%03d_truth.tsv.gz", args$anc1_name, sim_i)), delim = "\t")
    write_delim(pheno_2_data, file.path(truth_dir, sprintf("%s_sim%03d_truth.tsv.gz", args$anc2_name, sim_i)), delim = "\t")

    fam_1 <- append_phenotype_to_fam(sim_fam_1, pheno_1_data, sim_i)
    fam_2 <- append_phenotype_to_fam(sim_fam_2, pheno_2_data, sim_i)

    shared <- combined_effs$shared_eligible
    simulation_manifest[[sim_i]] <- tibble(
        simulation = sim_i,
        seed = current_seed,
        requested_h2 = herit,
        requested_rho = rho,
        n_effect_variants_anc1 = sum(combined_effs$eligible_1),
        n_effect_variants_anc2 = sum(combined_effs$eligible_2),
        n_shared_effect_variants = sum(shared),
        realized_effect_rho = cor(combined_effs$beta_std_1[shared], combined_effs$beta_std_2[shared]),
        realized_h2_anc1 = unique(pheno_1_data$realized_h2),
        realized_h2_anc2 = unique(pheno_2_data$realized_h2)
    )

    cat(sprintf("  [Sim %d] Phenotype appended to .fam file (Total columns: %d)\n", sim_i, ncol(fam_1)))
}

write_delim(bind_rows(simulation_manifest), file.path(truth_dir, "simulation_manifest.tsv"), delim = "\t")

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
