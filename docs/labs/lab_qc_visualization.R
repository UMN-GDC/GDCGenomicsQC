# Lab: Visualizing and Working with Basic QC Outputs
#
# Interact with output files from the Basic QC Pipeline using R and tidyverse.
# Load, inspect, visualize, and analyze QC metrics including missingness, MAF,
# heterozygosity, and sex checks.
#
# Prerequisites:
# - Completed the QC pipeline tutorial (QC outputs present in OUT_DIR/)
# - R (>= 4.0) with tidyverse installed: install.packages("tidyverse")
# - Set OUT_DIR to your pipeline output directory

conda_env <- '/scratch.global/R25_files/r25_r_env/lib/R/library'
.libPaths(unique(c(.Library.site, .Library)))
.libPaths(c(conda_env, .libPaths()))

library(tidyverse)

OUT_DIR <- "/scratch.global/GDC/r25outputs/toy"
subset <- "full"

# Load sample missingness (from initial QC)
smiss_path <- file.path(OUT_DIR, subset, "initial.smiss")
smiss <- read_tsv(smiss_path)
glimpse(smiss)
summary(smiss$F_MISS)

# Load variant missingness
vmiss_path <- file.path(OUT_DIR, subset, "initial.vmiss")
vmiss <- read_tsv(vmiss_path)
glimpse(vmiss)

# Load HWE and allele frequencies
hwe <- read_tsv(file.path(OUT_DIR, subset, "zoomhwe.hwe"))
glimpse(hwe)

freq <- read_tsv(file.path(OUT_DIR, subset, "MAF_check.afreq"))
glimpse(freq)

# 1. Sample missingness histogram
ggplot(smiss, aes(x = F_MISS * 100)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = c(2, 10), color = "red", linetype = "dashed") +
  labs(
    title = "Sample Missingness Distribution",
    x = "Missing Genotypes (%)",
    y = "Number of Samples"
  ) +
  theme_minimal()

# 2. Variant missingness histogram
ggplot(vmiss, aes(x = F_MISS * 100)) +
  geom_histogram(bins = 30, fill = "darkgreen", alpha = 0.7) +
  geom_vline(xintercept = 2, color = "red", linetype = "dashed") +
  labs(
    title = "Variant Missingness Distribution",
    x = "Missing Calls (%)",
    y = "Number of Variants"
  ) +
  theme_minimal()

# 3. Allele frequency distribution
freq %>%
  ggplot(aes(x = ALT_FREQS)) +
  geom_histogram() +
  theme_minimal()

freq %>%
  ggplot(aes(x = ALT_FREQS)) +
  geom_histogram() +
  theme_minimal() +
  facet_wrap(~`#CHROM`)

# 4. Identify samples with high missingness (>2% threshold)
high_missing <- smiss %>%
  filter(F_MISS > 0.02) %>%
  arrange(desc(F_MISS))
print(high_missing)
