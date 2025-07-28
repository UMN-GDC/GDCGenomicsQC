path_to_packages <- "/home/gdc/public/Ref/R"
.libPaths(c(path_to_packages, .libPaths()))  # Ensure R looks there for packages


install_if_missing <- function(pkg, lib = path_to_packages) {
  if (!requireNamespace(pkg, quietly = TRUE, lib.loc = lib)) {
    install.packages(pkg, lib = lib, dependencies = TRUE, repos = "https://cloud.r-project.org")
  }
  suppressMessages(library(pkg, character.only = TRUE, lib.loc = lib))
}

# Load libraries
install_if_missing("ggplot2")
install_if_missing("data.table")
install_if_missing("tidyverse")
# install_if_missing("dplyr")

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript plot_pca.R /path/to/WORK/PCA")
}

work_dir <- args[1]

# File paths
eigenvec_file <- file.path(work_dir, "PCA/merged_dataset_pca.eigenvec")
fam_file <- file.path(work_dir, "PCA/merged_common_bi.fam")
ancestry_file <- file.path(work_dir, "ancestry_estimation/data.txt")
save_dir <- file.path(work_dir, "PCA")

# Read ancestry labels
ancestry <- fread(ancestry_file, header = FALSE, col.names = c("FID", "IID", "Ancestry")) %>%
  mutate(across(c(FID, IID), as.character))

# Read PCA output
pca <- fread(eigenvec_file, header = FALSE)
colnames(pca) <- c("FID", "IID", paste0("PC", 1:(ncol(pca) - 2)))
pca <- pca %>%
  mutate(across(c(FID, IID), as.character))

# Read .fam file
fam <- fread(fam_file, header = FALSE)
colnames(fam) <- c("FID", "IID", "PID", "MID", "Sex", "Pheno")
fam <- fam[, .(FID, IID, Pheno)] %>%
  mutate(across(c(FID, IID), as.character))

# Read additional ancestry data
extra_ancestry_file <- "/home/gdc/shared/rfmix_ref/super_population_map_file.txt"
extra_ancestry <- fread(extra_ancestry_file, header = FALSE, col.names = c("IID", "Ancestry")) %>%
  mutate(IID = as.character(IID)) %>%
  mutate(FID = IID) %>%  # Create dummy FID to match structure
  select(FID, IID, Ancestry)


# Combine original and extra ancestry info
ancestry_all <- bind_rows(ancestry, extra_ancestry) %>%
  distinct(FID, IID, .keep_all = TRUE)  # Keep only one row per person

# Merge all info using the combined ancestry
data <- pca %>%
  left_join(ancestry_all, by = c("FID", "IID")) %>%
  left_join(fam, by = c("FID", "IID")) %>%
  mutate(Pheno = factor(Pheno, levels = c(1, 2), labels = c("Reference", "Study")))


# Function to plot pairs of PCs
plot_pc <- function(df, pcx, pcy) {
  ggplot(df, aes(x = .data[[pcx]], y = .data[[pcy]], color = Ancestry, shape = Pheno)) +
    geom_point(size = 2, alpha = 0.6) +
    theme_minimal() +
    labs(title = paste(pcx, "vs", pcy), x = pcx, y = pcy) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11)
    )
}


# Save PC1 vs PC2
png(file.path(save_dir, "PC1_vs_PC2.png"), width = 7, height = 6, units = "in", res = 300)
print(plot_pc(data, "PC1", "PC2"))
dev.off()

# Save PC1 vs PC3
png(file.path(save_dir, "PC1_vs_PC3.png"), width = 7, height = 6, units = "in", res = 300)
print(plot_pc(data, "PC1", "PC3"))
dev.off()

# Save PC2 vs PC3
png(file.path(save_dir, "PC2_vs_PC3.png"), width = 7, height = 6, units = "in", res = 300)
print(plot_pc(data, "PC2", "PC3"))
dev.off()
