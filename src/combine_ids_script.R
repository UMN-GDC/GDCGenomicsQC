library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]
name <- args[2]

# Read ancestry file, skipping headers
ancestry <- read.table("ancestry.bed", skip = 2, header = FALSE)
colnames(ancestry) <- c("sample", "AFR", "AMR", "EAS", "EUR", "SAS")

# Extract numeric part from ancestry ID (e.g., "ancestry101"  101)
ancestry$individual_index <- as.integer(gsub("ancestry", "", ancestry$sample))

# Read the index-to-ID mapping
map <- read.table("ancestry_index_map.tsv", header = TRUE)

# Merge by individual_index
merged <- merge(ancestry, map, by = "individual_index")

# Reorder columns: use individual_id instead of sample
final <- merged[, c("individual_id", "AFR", "AMR", "EAS", "EUR", "SAS")]
colnames(final) <- c("sample", "AFR", "AMR", "EAS", "EUR", "SAS")

# Save the result
write.table(final, paste0("ancestry_", name,".txt"), sep = "\t", quote = FALSE, row.names = FALSE)