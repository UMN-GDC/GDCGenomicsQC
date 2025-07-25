args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]          # directory path
name <- args[2]         # filename, e.g. "kinships.genome"

# ---- Load Libraries ----
if (!requireNamespace("GWASTools", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org/")
  BiocManager::install("GWASTools")
}
suppressPackageStartupMessages({
  library(GWASTools)
  library(reshape2)
})

# ---- Read PLINK .genome File ----
kin <- read.table(file.path(dir, name), header = TRUE)
cat("Loaded kin data with", nrow(kin), "rows and", ncol(kin), "columns\n")
print(colnames(kin))

# ---- IBD Plot with GWASTools ----
pdf(file.path(dir, "ibd_plot.pdf"), width = 8, height = 6)
ibdPlot(
  k0 = kin$Z0,
  k1 = kin$Z1,
  k2 = kin$Z2,
  kinship = kin$PI_HAT / 2
)
dev.off()
cat("IBD plot saved to ibd_plot.pdf\n")

# ---- Kinship Matrix Construction ----
kin_reformat <- data.frame(
  ID1 = paste(kin$FID1, kin$IID1),
  ID2 = paste(kin$FID2, kin$IID2),
  Kinship = kin$PI_HAT / 2
)

kin_matrix <- dcast(kin_reformat, ID1 ~ ID2, value.var = "Kinship")
rownames(kin_matrix) <- kin_matrix$ID1
kin_matrix <- as.matrix(kin_matrix[, -1])

# Save kinship matrix to RDS
saveRDS(kin_matrix, file = file.path(dir, "kinship_matrix.rds"))
cat("Kinship matrix saved to kinship_matrix.rds\n")

# ---- Identify Related Individuals (Kinship > 0.05) ----
coords <- which(kin_matrix > 0.05, arr.ind = TRUE)
row_ids <- rownames(kin_matrix)[coords[, 1]]
col_ids <- colnames(kin_matrix)[coords[, 2]]
ids_all <- c(row_ids, col_ids)
ids_unique <- unique(ids_all)

# ---- Output IDs to Exclude ----
exclude_file <- file.path(dir, "to_exclude.txt")
write.table(data.frame(ids_unique, ids_unique), file = exclude_file,
            col.names = FALSE, row.names = FALSE, quote = FALSE)

cat("List of individuals to exclude saved to to_exclude.txt\n")