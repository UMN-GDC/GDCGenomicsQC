# Load required libraries
library(dplyr)
library(readr)

# Define input directory and output directory
args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]   # Change this to your directory
output_dir <- args[2]   # Change this to where you want to save BED files 

# Ensure output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Get chromosome files (modify pattern if needed)
chromosomes <- 1:22  # Adjust if needed

# Loop through each chromosome
for (chr in chromosomes) {
  # Define file path
  msp_file <- file.path(input_dir, paste0("ancestry_chr", chr, ".msp.tsv"))
  sis_file <- file.path(input_dir, paste0("ancestry_chr", chr, ".sis.tsv"))

  # Skip if file does not exist
  if (!file.exists(msp_file) | !file.exists(sis_file)) {
    next
  }

  # Read Sample IDs from .sis.tsv
  sample_names <- read_delim(sis_file, delim = "\t", col_names = FALSE)[[1]]

  # Read local ancestry data
  msp_data <- read_delim(msp_file, delim = "\t", skip = 1, col_names = FALSE)

  # Assign column names
  colnames(msp_data) <- c("Chromosome", "Start", "End", "Marker", "Dummy", sample_names)

  # Remove the dummy column (column 5)
  msp_data <- msp_data %>% select(-Dummy)

  # Extract and save ancestry data for each individual
  for (sample in sample_names) {
    # Select relevant columns
    sample_data <- msp_data %>% select(Chromosome, Start, End, sample)

    # Define output file path
    output_bed <- file.path(output_dir, paste0(sample, "_chr", chr, ".bed"))

    # Save as BED file (no header)
    write_delim(sample_data, output_bed, delim = "\t", col_names = FALSE)
    
    print(paste("Saved:", output_bed))
  }
}
