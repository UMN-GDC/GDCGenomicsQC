#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
WORK <- args[1]
NAME <- args[2]

library(SNPRelate)
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(GENESIS)

# ---- Step 1: File paths ----
plink_prefix <- file.path(WORK, "Initial_QC", "QC4")
gds_file <- file.path(WORK, "gds", paste0(NAME, ".gds"))
seq_gds_file <- sub("\\.gds$", "_seq.gds", gds_file)
out_dir <- file.path(WORK, "pca_ir")
dir.create(dirname(gds_file), showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Step 2: Convert PLINK to SNPRelate GDS ----
snpgdsBED2GDS(
  bed.fn = paste0(plink_prefix, ".bed"),
  bim.fn = paste0(plink_prefix, ".bim"),
  fam.fn = paste0(plink_prefix, ".fam"),
  out.gdsfn = gds_file
)

# ---- Step 3: Open GDS file ----
genofile <- snpgdsOpen(gds_file)

# ---- Step 4: Run KING-robust relatedness ----
kingmat <- snpgdsIBDKING(genofile, verbose = TRUE)
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
kingmat2 <- kingToMatrix(kingmat, sample.id)

# ---- Step 5: PC-AiR ----
pcaobj <- pcair(genofile, kinobj = kingmat2, divobj = NULL, verbose = TRUE)

# ---- Step 6: Convert to SeqArray GDS for PC-Relate ----
# Convert SNPRelate GDS to SeqArray format
snpgdsClose(genofile)
seqSNP2GDS(gds.fn = gds_file, out.fn = seq_gds_file, verbose = TRUE)

# Open the resulting SeqArray file
seqfile <- seqOpen(seq_gds_file)
seqData <- SeqVarData(seqfile)

# ---- Step 7: PC-Relate ----
## -- how many variants?  -----------------------------------------------------
## 1) pull the variant IDs from the file …
# keep every variant
variant.id     <- seqGetData(seqfile, "variant.id")
variantFilter  <- list(variant.sel = rep(TRUE, length(variant.id)))

# now build the iterator
seqIter <- SeqVarIterator(seqData, variantFilter = variantFilter)


## -- PC-Relate ---------------------------------------------------------------
relate <- pcrelate(
  seqIter,
  pcs           = pcaobj$vectors,
  training.set  = pcaobj$unrels,
  verbose       = TRUE
)

# ---- Step 8: Save results ----
saveRDS(pcaobj, file = file.path(out_dir, paste0(NAME, "_pcaobj.RDS")))
saveRDS(relate, file = file.path(out_dir, paste0(NAME, "_pcrelate.RDS")))
write.table(pcaobj$unrels,
            file = file.path(out_dir, paste0(NAME, "_unrelated_ids.txt")),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# ---- Step 9: Cleanup ----
seqClose(seqfile)
