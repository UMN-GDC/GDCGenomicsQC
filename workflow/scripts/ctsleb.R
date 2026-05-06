library(CTSLEB)
library(data.table)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

# Parse arguments
sumstats_file <- NA
bed_file <- NA
bim_file <- NA
fam_file <- NA
pheno_file <- NA
out_dir <- "."

i <- 1
while (i <= length(args)) {
  if (args[i] == "--sumstats") {
    sumstats_file <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--bed") {
    bed_file <- sub("\\.bed$", "", args[i + 1])
    i <- i + 2
  } else if (args[i] == "--bim") {
    bim_file <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--fam") {
    fam_file <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--pheno") {
    pheno_file <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--out-dir") {
    out_dir <- args[i + 1]
    i <- i + 2
  } else {
    i <- i + 1
  }
}

if (is.na(sumstats_file) || is.na(bed_file)) {
  stop("Usage: Rscript ctsleb.R --sumstats <file> --bed <plink_prefix> --bim <file> --fam <file> --pheno <file> --out-dir <dir>")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Load summary statistics
sumstats <- fread(sumstats_file, header = TRUE)

# Set up PRS farm parameters
plink19_exec <- "plink"
plink2_exec <- "plink2"

PRS_farm <- SetParamsFarm(
  plink19_exec = plink19_exec,
  plink2_exec = plink2_exec,
  mem = 60000
)

# Run CTSLEB
prs_mat <- dimCT(
  results_dir = out_dir,
  sum_target = sumstats,
  sum_ref = sumstats,
  ref_plink = bed_file,
  target_plink = bed_file,
  test_target_plink = bed_file,
  out_prefix = "ctsleb",
  params_farm = PRS_farm
)

print("CTSLEB completed successfully")
