args <- commandArgs(trailingOnly = TRUE)
dir  <- args[1]                 # directory path
name <- args[2]                 # KING prefix, e.g. "king" (expects king.kin0, king.kin)

# ---- Libraries ----
if (!requireNamespace("GWASTools", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org/")
  BiocManager::install("GWASTools")
}
suppressPackageStartupMessages({
  library(GWASTools)
  library(reshape2)
})

# ---- Helper: read a KING table safely ----
read_king <- function(path) {
  if (!file.exists(path)) return(NULL)
  df <- tryCatch(read.table(path, header = TRUE, stringsAsFactors = FALSE),
                 error = function(e) NULL)
  if (is.null(df)) return(NULL)
  
  # Normalize column names across KING versions
  # Expected: FID1, ID1 (or IID1), FID2, ID2 (or IID2), NSNP/N_SNP, HetHet/HETHET, IBS0, Kinship
  cn <- colnames(df)
  if ("IID1" %in% cn) names(df)[names(df)=="IID1"] <- "ID1"
  if ("IID2" %in% cn) names(df)[names(df)=="IID2"] <- "ID2"
  if ("N_SNP" %in% cn) names(df)[names(df)=="N_SNP"] <- "NSNP"
  if ("HETHET" %in% cn) names(df)[names(df)=="HETHET"] <- "HetHet"
  if ("KINSHIP" %in% cn) names(df)[names(df)=="KINSHIP"] <- "Kinship"
  if ("IBS0" %in% cn) names(df)[names(df)=="IBS0"] <- "IBS0"
  
  # Keep only needed columns if present
  keep <- intersect(c("FID1","ID1","FID2","ID2","NSNP","HetHet","IBS0","Kinship"), colnames(df))
  df <- df[, keep, drop = FALSE]
  
  # Coerce numerics
  for (nm in c("NSNP","HetHet","IBS0","Kinship")) {
    if (nm %in% names(df)) df[[nm]] <- suppressWarnings(as.numeric(df[[nm]]))
  }
  df
}

# ---- Read KING outputs (kin0 + kin if present) ----
kin0_path <- file.path(dir, paste0(name, ".kin0"))
kin_path  <- file.path(dir, paste0(name, ".kin"))

kin0 <- read_king(kin0_path)
kinw <- read_king(kin_path)

if (is.null(kin0) && is.null(kinw)) {
  stop("No KING files found: ", kin0_path, " or ", kin_path)
}

kin <- do.call(rbind, Filter(Negate(is.null), list(kin0, kinw)))
cat("Loaded KING kinship pairs:", nrow(kin), "rows\n")
print(colnames(kin))

# ---- IBD-like plot using KING metrics ----
# GWASTools::ibdPlot can take KING kinship (kc) and IBS0 directly.
# (kc = KING kinship coefficient; ibs0 = opposite-homozygote proportion)
# Ref: GWASTools manual ibdPlot args. 
pdf(file.path(dir, "ibd_plot_king.pdf"), width = 8, height = 6)
ibdPlot(
  kc   = kin$Kinship,
  ibs0 = kin$IBS0,
  alpha = 0.05
)
dev.off()
cat("IBD/KING plot saved to ibd_plot_king.pdf\n")  # :contentReference[oaicite:2]{index=2}

# ---- Optional: assign degree categories from KING kc + IBS0 ----
# Default KING cutpoints (GWASTools) for kc:
# dup >= 0.3536; FS/PO >= 0.1768; deg2 >= 0.0884; deg3 >= 0.0442
# We'll attach a simple label for convenience.
assign_deg <- function(kc) {
  ifelse(kc >= 1/(2^(3/2)),  "dup/MZ",
         ifelse(kc >= 1/(2^(5/2)),  "1st",
                ifelse(kc >= 1/(2^(7/2)),  "2nd",
                       ifelse(kc >= 1/(2^(9/2)),  "3rd", "unrelated"))))
}
kin$degree <- assign_deg(kin$Kinship)
cat("Pair counts by category:\n")
print(table(kin$degree))  # :contentReference[oaicite:3]{index=3}

# ---- Kinship matrix (KING kinship coefficient) ----
id1 <- if ("FID1" %in% names(kin)) paste(kin$FID1, kin$ID1) else kin$ID1
id2 <- if ("FID2" %in% names(kin)) paste(kin$FID2, kin$ID2) else kin$ID2

kin_reformat <- data.frame(
  ID1 = id1,
  ID2 = id2,
  Kinship = kin$Kinship
)

kin_matrix <- dcast(kin_reformat, ID1 ~ ID2, value.var = "Kinship")
rownames(kin_matrix) <- kin_matrix$ID1
kin_matrix <- as.matrix(kin_matrix[, -1, drop = FALSE])
# Symmetrize: ensure missing cells are 0 on off-diagonals, 1/2 on diagonal if desired
# (KING kinship for identical samples is 0.5; we leave diagonal as 0/NA since pairs file lacks self-pairs)
kin_matrix[is.na(kin_matrix)] <- 0

saveRDS(kin_matrix, file = file.path(dir, "kinship_matrix.rds"))
cat("Kinship matrix saved to kinship_matrix.rds\n")

# ---- Identify related individuals (default: >= 3rd-degree threshold) ----
kc_cutoff <- 1/(2^(9/2))   # 0.04419417...  (change to 0.0884 for 2nd-degree)
rel_pairs <- subset(kin_reformat, Kinship >= kc_cutoff)

# Unique IDs to exclude (simple union of any pair above threshold)
ids_unique <- sort(unique(c(rel_pairs$ID1, rel_pairs$ID2)))

# Write two flavors:
# 1) Proper two-column FID IID (if FIDs present)
if ("FID1" %in% names(kin)) {
  split_ids <- do.call(rbind, strsplit(ids_unique, " ", fixed = TRUE))
  exclude_fid_iid <- data.frame(FID = split_ids[,1], IID = split_ids[,2], stringsAsFactors = FALSE)
  write.table(exclude_fid_iid,
              file = file.path(dir, "to_exclude.fid_iid.txt"),
              col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
  cat("Exclude list (FID IID) saved to to_exclude.fid_iid.txt\n")
}

# 2) Backward-compatible two identical columns of the combined ID (as in your original script)
write.table(data.frame(ids_unique, ids_unique),
            file = file.path(dir, "to_exclude.txt"),
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
cat("Exclude list (duplicated combined IDs) saved to to_exclude.txt\n")

cat(sprintf("Using cutoff kc >= %.4f yielded %d unique IDs to exclude.\n",
            kc_cutoff, length(ids_unique)))

# args <- commandArgs(trailingOnly = TRUE)
# dir <- args[1]          # directory path
# name <- args[2]         # filename, e.g. "kinships.genome"
# 
# # ---- Load Libraries ----
# if (!requireNamespace("GWASTools", quietly = TRUE)) {
#   install.packages("BiocManager", repos = "https://cloud.r-project.org/")
#   BiocManager::install("GWASTools")
# }
# suppressPackageStartupMessages({
#   library(GWASTools)
#   library(reshape2)
# })
# 
# # ---- Read PLINK .genome File ----
# kin <- read.table(file.path(dir, name), header = TRUE)
# cat("Loaded kin data with", nrow(kin), "rows and", ncol(kin), "columns\n")
# print(colnames(kin))
# 
# # ---- IBD Plot with GWASTools ----
# pdf(file.path(dir, "ibd_plot.pdf"), width = 8, height = 6)
# ibdPlot(
#   k0 = kin$Z0,
#   k1 = kin$Z1,
#   k2 = kin$Z2,
#   kinship = kin$PI_HAT / 2
# )
# dev.off()
# cat("IBD plot saved to ibd_plot.pdf\n")
# 
# # ---- Kinship Matrix Construction ----
# kin_reformat <- data.frame(
#   ID1 = paste(kin$FID1, kin$IID1),
#   ID2 = paste(kin$FID2, kin$IID2),
#   Kinship = kin$PI_HAT / 2
# )
# 
# kin_matrix <- dcast(kin_reformat, ID1 ~ ID2, value.var = "Kinship")
# rownames(kin_matrix) <- kin_matrix$ID1
# kin_matrix <- as.matrix(kin_matrix[, -1])
# 
# # Save kinship matrix to RDS
# saveRDS(kin_matrix, file = file.path(dir, "kinship_matrix.rds"))
# cat("Kinship matrix saved to kinship_matrix.rds\n")
# 
# # ---- Identify Related Individuals (Kinship > 0.05) ----
# coords <- which(kin_matrix > 0.05, arr.ind = TRUE)
# row_ids <- rownames(kin_matrix)[coords[, 1]]
# col_ids <- colnames(kin_matrix)[coords[, 2]]
# ids_all <- c(row_ids, col_ids)
# ids_unique <- unique(ids_all)
# 
# # ---- Output IDs to Exclude ----
# exclude_file <- file.path(dir, "to_exclude.txt")
# write.table(data.frame(ids_unique, ids_unique), file = exclude_file,
#             col.names = FALSE, row.names = FALSE, quote = FALSE)
# 
# cat("List of individuals to exclude saved to to_exclude.txt\n")
