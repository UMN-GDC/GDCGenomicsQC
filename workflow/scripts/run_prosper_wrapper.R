#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
})

parser <- ArgumentParser(description = "Run PROSPER multi-ancestry PRS workflow")
parser$add_argument("--package", required = TRUE, help = "Path to PROSPER package directory")
parser$add_argument("--out-dir", required = TRUE, help = "Output directory")
parser$add_argument("--plink2", required = TRUE, help = "Path to plink2 executable")
parser$add_argument("--target-sumstats", required = TRUE)
parser$add_argument("--aux-sumstats", required = TRUE)
parser$add_argument("--pop", required = TRUE, help = "Comma-separated population labels in the same order as summary stats")
parser$add_argument("--lassosum-param-files", required = TRUE, help = "Comma-separated lassosum2 tuning parameter files in the same order as summary stats")
parser$add_argument("--prefix", required = TRUE, help = "Target-population output prefix")
parser$add_argument("--tuning-bfile", required = TRUE, help = "PLINK prefix for tuning set")
parser$add_argument("--testing-bfile", required = TRUE, help = "PLINK prefix for testing set")
parser$add_argument("--pheno-tuning", default = NULL)
parser$add_argument("--covar-tuning", default = NULL)
parser$add_argument("--pheno-testing", default = NULL)
parser$add_argument("--covar-testing", default = NULL)
parser$add_argument("--sl-library", default = "SL.glmnet,SL.ridge,SL.lm")
parser$add_argument("--linear-score", default = "TRUE")
parser$add_argument("--cleanup", default = "TRUE")
parser$add_argument("--chrom", default = "1-22")
parser$add_argument("--ll", default = "5")
parser$add_argument("--lc", default = "5")
parser$add_argument("--ncores", default = "1")
parser$add_argument("--verbose", default = "1")
args <- parser$parse_args()

dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)

prosper_script <- file.path(args$package, "scripts", "PROSPER.R")
tuning_script <- file.path(args$package, "scripts", "tuning_testing.R")

if (!file.exists(prosper_script)) {
  stop("Could not find PROSPER.R at: ", prosper_script)
}
if (!file.exists(tuning_script)) {
  stop("Could not find tuning_testing.R at: ", tuning_script)
}

required_files <- c(
  args$target_sumstats,
  args$aux_sumstats,
  strsplit(args$lassosum_param_files, ",", fixed = TRUE)[[1]]
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    "Missing required PROSPER input files:\n",
    paste(missing_files, collapse = "\n"),
    "\nGenerate the lassosum2 parameter files first and pass them via --lassosum-param-files."
  )
}


run_rscript <- function(script, args_vec) {
  cmd <- c(script, args_vec)
  rscript <- file.path(R.home("bin"), "Rscript")

  prosper_bin <- "/scratch.global/saonli/conda-envs/prosper/bin"
  Sys.setenv(PATH = paste(prosper_bin, Sys.getenv("PATH"), sep = ":"))

  makevars_dir <- tempfile("r_makevars_")
  dir.create(makevars_dir, recursive = TRUE, showWarnings = FALSE)
  makevars <- file.path(makevars_dir, "Makevars")
  writeLines(c(
    "CC=gcc",
    "CXX=g++",
    "CXX11=g++",
    "CXX14=g++",
    "CXX17=g++",
    "CXX20=g++"
  ), makevars)
  Sys.setenv(R_MAKEVARS_USER = makevars)

  message("Running: ", rscript, " ", paste(cmd, collapse = " "))
  status <- system2(rscript, cmd)
  if (!identical(status, 0L)) {
    stop("Failed while running ", basename(script), " (exit code ", status, ")")
  }
}



prosper_args <- c(
  "--PATH_package", args$package,
  "--PATH_out", args$out_dir,
  "--FILE_sst", paste(args$target_sumstats, args$aux_sumstats, sep = ","),
  "--pop", args$pop,
  "--lassosum_param", args$lassosum_param_files,
  "--chrom", args$chrom,
  "--Ll", args$ll,
  "--Lc", args$lc,
  "--verbose", args$verbose,
  "--NCORES", args$ncores
)

tuning_args <- c(
  "--PATH_out", args$out_dir,
  "--PATH_plink", args$plink2,
  "--prefix", args$prefix,
  "--SL_library", args$sl_library,
  "--linear_score", args$linear_score,
  "--bfile_tuning", args$tuning_bfile,
  "--testing", "TRUE",
  "--bfile_testing", args$testing_bfile,
  "--verbose", args$verbose,
  "--cleanup", args$cleanup,
  "--NCORES", args$ncores
)

if (!is.null(args$pheno_tuning) && nzchar(args$pheno_tuning)) {
  tuning_args <- c(tuning_args, "--pheno_tuning", args$pheno_tuning)
}
if (!is.null(args$covar_tuning) && nzchar(args$covar_tuning)) {
  tuning_args <- c(tuning_args, "--covar_tuning", args$covar_tuning)
}
if (!is.null(args$pheno_testing) && nzchar(args$pheno_testing)) {
  tuning_args <- c(tuning_args, "--pheno_testing", args$pheno_testing)
}
if (!is.null(args$covar_testing) && nzchar(args$covar_testing)) {
  tuning_args <- c(tuning_args, "--covar_testing", args$covar_testing)
}

run_rscript(prosper_script, prosper_args)
run_rscript(tuning_script, tuning_args)

after_dir <- file.path(args$out_dir, paste0("after_ensemble_", args$prefix))
r2_file <- file.path(after_dir, "R2.txt")
prs_file <- file.path(after_dir, "PROSPER_prs_file.txt")

if (!file.exists(r2_file)) {
  stop("PROSPER finished but did not create expected R2 file: ", r2_file)
}

if (tolower(args$linear_score) == "true" && !file.exists(prs_file)) {
  stop("PROSPER finished but did not create expected PRS file: ", prs_file)
}

message("PROSPER completed successfully.")
message("Final PRS file: ", prs_file)
message("Performance summary: ", r2_file)
