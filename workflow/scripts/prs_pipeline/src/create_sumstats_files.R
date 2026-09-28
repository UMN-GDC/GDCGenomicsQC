#!/usr/bin/env Rscript

# -------------------------------------------------------------------------
# Setup
# -------------------------------------------------------------------------
.libPaths("/projects/standard/gdc/public/Ref/R")
suppressMessages(library(readr))
suppressMessages(library(tidyverse))

# -------------------------------------------------------------------------
# Helper Functions
# -------------------------------------------------------------------------

# Read association file
read_assoc <- function(path) {
  read_table(path, col_types = cols()) %>%
    transmute(
      rsid = SNP,
      chr = CHR,
      effect_allele_assoc = A1,
      beta = BETA,
      beta_se = SE,
      p = P,
      n_eff = n_eff
    )
}

# Read PLINK .bim file
read_bim <- function(path) {
  read_table(path, col_names = FALSE, col_types = cols()) %>%
    set_names(c("chrom", "rsid", "cm", "bp", "a1_bim", "a0_bim")) %>%
    transmute(rsid, a1_bim, a0_bim, bp)
}

# Merge assoc + bim and flip betas where necessary
merge_with_bim <- function(assoc_df, bim_df, keep_a0 = TRUE) {
  assoc_df %>%
    inner_join(bim_df, by = "rsid") %>%
    mutate(
      flip = effect_allele_assoc != a1_bim,
      beta = ifelse(flip, -beta, beta),
      a1 = effect_allele_assoc,
      a0 = ifelse(flip, a1_bim, a0_bim)
    ) %>%
    {
      if (keep_a0) {
        select(., rsid, chr, a1, a0, beta, beta_se, n_eff)
      } else {
        select(., rsid, chr, a1, beta, beta_se, n_eff)
      }
    } %>%
    na.omit()
}

# Write CTSLEB format
write_ctsleb <- function(assoc_df, bim_df, out_path) {
  assoc_df %>%
    inner_join(bim_df, by = "rsid") %>%
    mutate(
      CHR = chr,
      SNP = rsid,
      BP = bp,
      A1 = effect_allele_assoc,
      rs_id = rsid
    ) %>%
    select(CHR, SNP, BP, A1, beta, beta_se, p, rs_id) %>%
    na.omit() %>%
    write.table(file = out_path, sep = " ", quote = FALSE,
                row.names = FALSE, col.names = TRUE)
}

# -------------------------------------------------------------------------
# Argument Parsing
# -------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (!(length(args) %in% c(4, 7))) {
  stop("Script requires 4 or 7 arguments.", call. = FALSE)
}

base_location          <- args[1]
simulated_assoc_name_1 <- args[2]
simulated_bim_name_1   <- args[3]
sum_stat_1_name        <- args[4]

has_training <- length(args) == 7
if (has_training) {
  simulated_assoc_name_2 <- args[5]
  simulated_bim_name_2   <- args[6]
  sum_stat_2_name        <- args[7]
}

setwd(base_location)

# -------------------------------------------------------------------------
# Target (Population 1) Output
# -------------------------------------------------------------------------

assoc1 <- read_assoc(simulated_assoc_name_1)
bim1   <- read_bim(simulated_bim_name_1)

# PROSPER-format (main sumstats)
prosper_sumstats <- merge_with_bim(assoc1, bim1, keep_a0 = TRUE)
write.table(prosper_sumstats, file = sum_stat_1_name,
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# CTSLEB version
target_ctsleb_name <- paste0(tools::file_path_sans_ext(sum_stat_1_name), "_CTSLEB.txt")
write_ctsleb(assoc1, bim1, target_ctsleb_name)

# -------------------------------------------------------------------------
# Training (Population 2) Output  optional
# -------------------------------------------------------------------------

if (has_training) {
  assoc2 <- read_assoc(simulated_assoc_name_2)
  bim2   <- read_bim(simulated_bim_name_2)

  # PROSPER-format
  training_sumstats <- merge_with_bim(assoc2, bim2, keep_a0 = TRUE)
  write.table(training_sumstats, file = sum_stat_2_name,
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

  # CTSLEB version
  training_ctsleb_name <- paste0(tools::file_path_sans_ext(sum_stat_2_name), "_CTSLEB.txt")
  write_ctsleb(assoc2, bim2, training_ctsleb_name)
}

