#!/usr/bin/env Rscript

library(argparse)
library(dplyr)
library(bigreadr)

# 1. SET UP ARGUMENT PARSER
parser <- ArgumentParser(description='Clean GWAS summary stats and calculate n_eff')
parser$add_argument("--input", type="character", required=TRUE, help="Path to tab-separated summary stats")
parser$add_argument("--bim", type="character", required=TRUE, help="Path to genotype .bim file for allele alignment")
parser$add_argument("--output", type="character", required=TRUE, help="Path for the output file")
parser$add_argument("--n_total", type="integer", help="Total sample size (if constant)")

args <- parser$parse_args()

# 2. LOAD DATA and ensure we have proper inputs
message("Reading input file: ", args$input)
df <- bigreadr::fread2(args$input)

target_rsid <- c("rsid", "rs_id", "rsids")
actual_col <- intersect(target_rsid, names(df))[1]
if (is.na(actual_col)) {
  stop("Error: Could not find an RSID column. Expected one of: ", 
       paste(target_rsid, collapse = ", "))
}

a1_options <- c("A1", "alt", "a1")
a1_present <- intersect(a1_options, names(df))[1]
if (is.na(a1_present)) {
  stop("Error: Could not find an a1 column. Expected one of: ", 
       paste(a1_options, collapse = ", "))
}

pval_options <- c("p", "pval")
pval_present <- intersect(pval_options, names(df))[1]
if (is.na(pval_present)) {
  stop("Error: Could not find a pval column. Expected one of: ", 
       paste(pval_options, collapse = ", "))
}

sebeta_options <- c("sebeta", "beta_se")
sebeta_present <- intersect(sebeta_options, names(df))[1]
if (is.na(sebeta_present)) {
  stop("Error: Could not find a sebeta column. Expected one of: ", 
       paste(sebeta_options, collapse = ", "))
}

message("Initial rows: ", nrow(df))

df <- df %>%
  rename(
    rsid = all_of(actual_col),
    alt = all_of(a1_present),
    pval = all_of(pval_present),
    sebeta = all_of(sebeta_present)
  )

message("Reading BIM file: ", args$bim)
bim <- bigreadr::fread2(args$bim, select = c(1, 2, 4, 5, 6))
colnames(bim) <- c("bim.chr", "rsid", "bim.pos", "bim.a1", "bim.a0")

# Cleaning files
message("Cleaning RSID columns for matching...")
df <- df %>%
  mutate(
    rsid = sub(",.*", "", rsid),
    rsid = trimws(rsid)
  )

bim <- bim %>%
  mutate(rsid = trimws(rsid))

message("Aligning alleles by RSID and processing...")

processed_df <- df %>%
  # Join with BIM to ensure genomic coordinates match your genotypes
  inner_join(bim, by = c("rsid")) %>%
  mutate(
    n_eff = if (!is.null(args$n_total)) {
      args$n_total
    } else {
      1 / (2 * af * (1 - af) * sebeta^2)
    },

    # Standardize Reference Allele (a0)
    a1 = alt,
    a0 = ifelse(a1 == bim.a1, bim.a0, bim.a1),

    # Ensure P-value is present
    CHR = as.integer(bim.chr),
    BP = as.integer(bim.pos),
    P = pval
  ) %>%
  # Format specifically for PRS tools
  select(
    SNP = rsid,
    CHR = CHR,
    BP = BP,
    A1 = a1,
    A2 = a0,
    beta = beta,
    beta_se = sebeta,
    P = P,
    n_eff = n_eff
  ) %>%
  # Remove any rows where n_eff calculation resulted in Inf or NA due to AF=0 or 1
  filter(is.finite(n_eff)) %>%
  # CRITICAL FIX: Drop duplicate variants so PLINK --score doesn't crash
  distinct(SNP, .keep_all = TRUE)

message("Rows after cleaning, joining, and deduplication: ", nrow(processed_df))

# 4. SAVE OUTPUT
message("Writing cleaned stats to: ", args$output)
write.table(processed_df, 
            file = args$output, 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)

message("Done! Processed ", nrow(processed_df), " unique SNPs.")