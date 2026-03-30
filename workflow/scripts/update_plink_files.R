#!/usr/bin/env Rscript

# ==========================================================
# Title: Generate Updated PLINK Files for Simulated Populations
# Author: Kody DeGolier
# Description:
#   This script prepares simulated phenotype data for two ancestries,
#   applies common SNP filters, and generates updated PLINK datasets.
# Usage:
#   Rscript script_name.R <location> <anc1_sim_file> <anc1_plink> 
#                         <anc2_sim_file> <anc2_plink> 
#                         <anc1_name> <anc2_name> <common_snps_file>
# Example:
#   Rscript prepare_plink.R ~/data AFR_sim.tsv 1kg_AFR EUR_sim.tsv 1kg_EUR AFR EUR common_snps.txt
# ==========================================================

suppressMessages(library(tidyverse))

# -------------------------------
# Parse Command-Line Arguments
# -------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 8) {
  stop("Usage: Rscript script.R <location> <anc1_sim> <anc1_plink> <anc2_sim> <anc2_plink> <anc1> <anc2> <common_snps>", call. = FALSE)
}

location            <- args[1]
anc1_simulated_file <- args[2]
anc1_plink_file     <- args[3]
anc2_simulated_file <- args[4]
anc2_plink_file     <- args[5]
anc1                <- args[6]
anc2                <- args[7]
common_snps_file    <- args[8]

# -------------------------------
# Basic Input Validation
# -------------------------------
required_files <- c(anc1_simulated_file, anc2_simulated_file, common_snps_file)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(paste("The following files were not found:\n", paste(missing_files, collapse = "\n")), call. = FALSE)
}

if (!dir.exists(location)) {
  stop(paste("The provided directory does not exist:", location), call. = FALSE)
}

# -------------------------------
# Helper Function for Logging
# -------------------------------
log_msg <- function(...) cat("[INFO]", ..., "\n")

# -------------------------------
# Initialize and Print Summary
# -------------------------------
setwd(location)
log_msg("Working directory set to:", location)
log_msg("Ancestry 1:", anc1)
log_msg("Ancestry 2:", anc2)
log_msg("Common SNP list:", common_snps_file)

# -------------------------------
# Step 1: Read and Update Simulated Data
# -------------------------------
process_sim_file <- function(sim_file, suffix) {
  log_msg("Reading simulated phenotype file:", sim_file)
  data <- read.table(sim_file, header = TRUE, sep = "\t", check.names = FALSE)
  # data$FID <- 0
  file_root <- sub("\\.[^\\.]*$", "", sim_file)
  output_file <- paste0(file_root, "_updated_", suffix, ".txt")
  log_msg("Writing updated phenotype file:", output_file)
  write.table(data, file = output_file, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  return(output_file)
}

anc1_updated_sim_file <- process_sim_file(anc1_simulated_file, "1")
anc2_updated_sim_file <- process_sim_file(anc2_simulated_file, "2")

# -------------------------------
# Step 2: Run PLINK Commands
# -------------------------------
run_plink <- function(plink_base, pheno_file, out_prefix) {
  cmd <- paste0(
    "plink2 --pfile ", plink_base,
    " --pheno ", pheno_file,
    " --pheno-name ", "pheno",
    " --make-bed --out ", out_prefix
  )
  log_msg("Running PLINK command:\n", cmd)
  status <- system(cmd, ignore.stdout = FALSE, ignore.stderr = FALSE)
  if (status != 0) stop(paste("PLINK command failed for", out_prefix))
}

final_name_anc1 <- file.path(location, paste0(anc1, "_simulation"))
final_name_anc2 <- file.path(location, paste0(anc2, "_simulation"))

run_plink(anc1_plink_file, anc1_simulated_file, final_name_anc1)
run_plink(anc2_plink_file, anc2_simulated_file, final_name_anc2)

# -------------------------------
# Step 3: Verify .fam Files
# -------------------------------
check_fam <- function(fam_path) {
  if (!file.exists(fam_path)) stop(paste("Missing .fam file:", fam_path))
  fam <- read.table(fam_path, header = FALSE, col.names = c("FID", "IID", "V3", "V4", "V5", "Phenotype"))
  log_msg("Preview of .fam file:", fam_path)
  print(head(fam))
  invisible(fam)
}

log_msg("Verifying output .fam files...")
check_fam(paste0(final_name_anc1, ".fam"))
check_fam(paste0(final_name_anc2, ".fam"))

# -------------------------------
# Completion
# -------------------------------
log_msg("All steps completed successfully.")
