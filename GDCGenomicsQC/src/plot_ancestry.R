# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

# Load RFMIX2 ancestry output
rfmix_output <- "local_ancestry_output.rfmix2.local_ancestry.txt"

# Read the data
ancestry_data <- read_delim(rfmix_output, delim = "\t", col_names = FALSE)

# Assign column names (modify based on your data structure)
colnames(ancestry_data) <- c("Sample", "Chromosome", "Position", "Ancestry")

# Count occurrences of each ancestry per chromosome
ancestry_counts <- ancestry_data %>%
  group_by(Chromosome, Ancestry) %>%
  summarise(Count = n(), .groups = "drop") %>%
  spread(Ancestry, Count, fill = 0)  # Convert to wide format

# Normalize for proportion
ancestry_counts <- ancestry_counts %>%
  mutate(across(-Chromosome, ~ . / sum(.), .names = "prop_{col}")) %>%
  pivot_longer(cols = starts_with("prop_"), names_to = "Ancestry", values_to = "Proportion") %>%
  mutate(Ancestry = gsub("prop_", "", Ancestry))  # Clean up ancestry labels

# Plot ancestry proportions per chromosome
ggplot(ancestry_counts, aes(x = as.factor(Chromosome), y = Proportion, fill = Ancestry)) +
  geom_bar(stat = "identity") +
  labs(title = "Ancestry Proportion per Chromosome",
       x = "Chromosome",
       y = "Proportion") +
  theme_minimal() +
  scale_fill_viridis_d()  # Use viridis color scheme

