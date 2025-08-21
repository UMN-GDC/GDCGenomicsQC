# use_king_related.R
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 2)
dir  <- args[1]
name <- args[2]
kc_cutoff <- if (length(args) >= 3) as.numeric(args[3]) else 1/(2^(9/2))  # 0.04419417...
ibs0_err  <- if (length(args) >= 4) as.numeric(args[4]) else 0.003        # KING PO IBS0 threshold

suppressPackageStartupMessages({
  library(reshape2)
})

# ---------- Robust reader for KING .kin0/.kin ----------
read_king <- function(path) {
  if (!file.exists(path)) return(NULL)
  df <- tryCatch(read.table(path, header = TRUE, stringsAsFactors = FALSE),
                 error = function(e) NULL)
  if (is.null(df)) return(NULL)
  cn <- colnames(df)
  if ("IID1" %in% cn) names(df)[names(df)=="IID1"] <- "ID1"
  if ("IID2" %in% cn) names(df)[names(df)=="IID2"] <- "ID2"
  if ("N_SNP" %in% cn) names(df)[names(df)=="N_SNP"] <- "NSNP"
  if ("HETHET" %in% cn) names(df)[names(df)=="HETHET"] <- "HetHet"
  if ("KINSHIP" %in% cn) names(df)[names(df)=="KINSHIP"] <- "Kinship"
  keep <- intersect(c("FID1","ID1","FID2","ID2","NSNP","HetHet","IBS0","Kinship"), names(df))
  df <- df[, keep, drop = FALSE]
  for (nm in intersect(c("NSNP","HetHet","IBS0","Kinship"), names(df))) {
    df[[nm]] <- suppressWarnings(as.numeric(df[[nm]]))
  }
  df
}

kin0 <- read_king(file.path(dir, paste0(name, ".kin0")))
kinw <- read_king(file.path(dir, paste0(name, ".kin")))
stopifnot(!is.null(kin0) || !is.null(kinw))
kin <- do.call(rbind, Filter(Negate(is.null), list(kin0, kinw)))

cat("Loaded KING pairs:", nrow(kin), "rows\n")
print(colnames(kin))

# ---------- KING classification (kc = Kinship, ibs0 = IBS0) ----------
# Thresholds from GWASTools docs: dup 0.3536, 1st 0.1768, 2nd 0.0884, 3rd 0.0442; PO separated by IBS0≈0. :contentReference[oaicite:2]{index=2}
cut_dup  <- 1/(2^(3/2))  # 0.3535534
cut_1st  <- 1/(2^(5/2))  # 0.1767767
cut_2nd  <- 1/(2^(7/2))  # 0.08838835
cut_3rd  <- 1/(2^(9/2))  # 0.04419417

kc <- kin$Kinship
ibs0 <- kin$IBS0

deg <- ifelse(kc >= cut_dup, "Dup/MZ",
       ifelse(kc >= cut_1st & ibs0 <= ibs0_err, "PO",
       ifelse(kc >= cut_1st, "FS",
       ifelse(kc >= cut_2nd, "2nd",
       ifelse(kc >= cut_3rd, "3rd", "U")))))
kin$degree <- deg
cat("Pair counts by KING category:\n"); print(table(kin$degree))

# ---------- Plot: IBS0 vs Kinship with KING cut lines ----------
pdf(file.path(dir, "king_scatter.pdf"), width = 8, height = 6)
plot(ibs0, kc, pch = 20, cex = 0.6, xlab = "IBS0 (opposite homozygote fraction)",
     ylab = "KING kinship coefficient", main = "KING: IBS0 vs kinship")
abline(h = cut_1st, lty = 2)  # 1st vs 2nd
abline(h = cut_2nd, lty = 2)  # 2nd vs 3rd
abline(h = cut_3rd, lty = 2)  # 3rd vs unrelated
abline(v = ibs0_err, lty = 3) # PO split (ibs0≈0)
legend("topright",
       legend = c(">=1st (0.1768)","2nd (0.0884)","3rd (0.0442)","PO split (IBS0)"),
       lty = c(2,2,2,3), bty = "n", cex = 0.9)
dev.off()
cat("Plot written: king_scatter.pdf\n")

# ---------- Kinship matrix ----------
id1 <- if ("FID1" %in% names(kin)) paste(kin$FID1, kin$ID1) else kin$ID1
id2 <- if ("FID2" %in% names(kin)) paste(kin$FID2, kin$ID2) else kin$ID2

kin_reformat <- data.frame(ID1 = id1, ID2 = id2, Kinship = kc)
kin_matrix <- dcast(kin_reformat, ID1 ~ ID2, value.var = "Kinship")
rownames(kin_matrix) <- kin_matrix$ID1
kin_matrix <- as.matrix(kin_matrix[, -1, drop = FALSE])
kin_matrix[is.na(kin_matrix)] <- 0
saveRDS(kin_matrix, file = file.path(dir, "kinship_matrix.rds"))
cat("Kinship matrix saved: kinship_matrix.rds\n")

# ---------- Exclude lists (kc >= kc_cutoff) ----------
rel_pairs <- subset(kin_reformat, Kinship >= kc_cutoff)
ids_unique <- sort(unique(c(rel_pairs$ID1, rel_pairs$ID2)))

# FID/IID two-column file (if FIDs present)
if ("FID1" %in% names(kin)) {
  split_ids <- do.call(rbind, strsplit(ids_unique, " ", fixed = TRUE))
  exclude_fid_iid <- data.frame(FID = split_ids[,1], IID = split_ids[,2], stringsAsFactors = FALSE)
  write.table(exclude_fid_iid, file = file.path(dir, "to_exclude.fid_iid.txt"),
              col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
  cat("Wrote: to_exclude.fid_iid.txt\n")
}

# Backward-compatible duplicate-ID format
write.table(data.frame(ids_unique, ids_unique),
            file = file.path(dir, "to_exclude.txt"),
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
cat(sprintf("kc_cutoff=%.5f ⇒ %d unique IDs to exclude. Wrote: to_exclude.txt\n",
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
