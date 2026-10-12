#!/usr/bin/env Rscript

# optparse (not argparse): this runs inside the upstream `singleprshelper`
# image, which ships optparse but not argparse.
suppressPackageStartupMessages({
  library(optparse)
})

opt_list <- list(
  make_option("--package", type = "character", help = "Path to PROSPER package directory"),
  make_option("--out-dir", type = "character", help = "Output directory"),
  make_option("--plink2", type = "character", default = "plink2", help = "Path to plink2 executable"),
  make_option("--target-sumstats", type = "character"),
  make_option("--aux-sumstats", type = "character"),
  make_option("--pop", type = "character", help = "Comma-separated population labels in the same order as summary stats"),
  make_option("--lassosum-param-files", type = "character", help = "Comma-separated lassosum2 tuning parameter files in the same order as summary stats"),
  make_option("--prefix", type = "character", help = "Target-population output prefix"),
  make_option("--tuning-bfile", type = "character", help = "PLINK prefix for tuning set"),
  make_option("--testing-bfile", type = "character", help = "PLINK prefix for testing set"),
  make_option("--pheno-tuning", type = "character", default = NULL),
  make_option("--covar-tuning", type = "character", default = NULL),
  make_option("--pheno-testing", type = "character", default = NULL),
  make_option("--covar-testing", type = "character", default = NULL),
  make_option("--sl-library", type = "character", default = "SL.glmnet,SL.ridge,SL.lm"),
  make_option("--linear-score", type = "logical", default = TRUE),
  make_option("--cleanup", type = "logical", default = TRUE),
  make_option("--chrom", type = "character", default = "1-22"),
  make_option("--ll", type = "character", default = "5"),
  make_option("--lc", type = "character", default = "5"),
  make_option("--ncores", type = "character", default = "1"),
  make_option("--verbose", type = "character", default = "1")
)
parser <- OptionParser(option_list = opt_list)
args <- parse_args(parser)
# This optparse build keeps dashes in dest names (e.g. "out-dir").
names(args) <- gsub("-", "_", names(args))

needed <- c("package", "out_dir", "target_sumstats", "aux_sumstats", "pop",
            "lassosum_param_files", "prefix", "tuning_bfile")
missing <- needed[unlist(lapply(needed, function(x) is.null(args[[x]])))]
if (length(missing) > 0) {
  stop("Missing required PROSPER options: ", paste(missing, collapse = ", "))
}

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

  # Optional compiler toolchain override (no hard-coded host paths): set
  # PROSPER_BIN to a directory to prepend to PATH, and PROSPER_MAKEVARS to a
  # Makevars path to point R at for sourceCpp. Inside singleprshelper both are
  # usually unset (the image ships its own compilers).
  prosper_bin <- Sys.getenv("PROSPER_BIN", "")
  if (nzchar(prosper_bin) && dir.exists(prosper_bin)) {
    Sys.setenv(PATH = paste(prosper_bin, Sys.getenv("PATH"), sep = ":"))
  }

  makevars <- Sys.getenv("PROSPER_MAKEVARS", "")
  if (nzchar(makevars)) {
    makevars_dir <- dirname(makevars)
    if (!dir.exists(makevars_dir)) {
      dir.create(makevars_dir, recursive = TRUE, showWarnings = FALSE)
    }
    writeLines(c(
      "CC=gcc",
      "CXX=g++",
      "CXX11=g++",
      "CXX14=g++",
      "CXX17=g++",
      "CXX20=g++"
    ), makevars)
    Sys.setenv(R_MAKEVARS_USER = makevars)
  }

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

if (isTRUE(args$linear_score) && !file.exists(prs_file)) {
  stop("PROSPER finished but did not create expected PRS file: ", prs_file)
}

message("PROSPER completed successfully.")
message("Final PRS file: ", prs_file)
message("Performance summary: ", r2_file)
