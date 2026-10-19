#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# Keep this wrapper runnable in the existing PROSPER R environment, which has
# data.table but not the optional argparse package. All CLI options take a
# value, so base R parsing is sufficient and avoids a separate dependency.
parse_arguments <- function(argv) {
  defaults <- list(
    pheno_tuning = NULL,
    covar_tuning = NULL,
    pheno_testing = NULL,
    covar_testing = NULL,
    sl_library = "SL.glmnet,SL.ridge,SL.lm",
    linear_score = "TRUE",
    cleanup = "TRUE",
    chrom = "1-22",
    ll = "5",
    lc = "5",
    ncores = "1",
    verbose = "1"
  )
  required <- c(
    "package", "out_dir", "plink2", "target_sumstats", "aux_sumstats",
    "pop", "lassosum_param_files", "prefix", "tuning_bfile", "testing_bfile"
  )
  values <- defaults
  index <- 1L

  while (index <= length(argv)) {
    option <- argv[[index]]
    if (!startsWith(option, "--") || index == length(argv)) {
      stop("Expected --option value pairs; received: ", option)
    }
    key <- gsub("-", "_", sub("^--", "", option))
    value <- argv[[index + 1L]]
    values[[key]] <- value
    index <- index + 2L
  }

  missing <- required[!vapply(required, function(key) {
    !is.null(values[[key]]) && nzchar(values[[key]])
  }, logical(1))]
  if (length(missing) > 0L) {
    stop("Missing required arguments: ", paste(paste0("--", gsub("_", "-", missing)), collapse = ", "))
  }
  values
}

args <- parse_arguments(commandArgs(trailingOnly = TRUE))

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

r2 <- fread(r2_file)
if (!all(c("tuning_R2", "testing_R2") %in% names(r2)) || nrow(r2) != 1L) {
  stop("Unexpected PROSPER R2.txt format; expected one row with tuning_R2 and testing_R2.")
}
testing_pheno <- if (!is.null(args$pheno_testing) && nzchar(args$pheno_testing)) {
  fread(args$pheno_testing)
} else {
  fread(paste0(args$testing_bfile, ".fam"))
}
fair_metrics <- data.table(
  method = "PROSPER",
  tuning_r2 = as.numeric(r2$tuning_R2[[1]]),
  test_r2 = as.numeric(r2$testing_R2[[1]]),
  test_n = nrow(testing_pheno)
)
fwrite(fair_metrics, file.path(after_dir, "PROSPER_fair_metrics.tsv"), sep = "\t")

message("PROSPER completed successfully.")
message("Final PRS file: ", prs_file)
message("Performance summary: ", r2_file)
