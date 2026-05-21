#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
    stop("Usage: Rscript run_pcair_pcrelate.R <bed_prefix> <out_dir> <gds_file> <seq_gds_file>\n")
}

bed_prefix <- args[1]
out_dir <- args[2]
gds_file <- args[3]
seq_gds_file <- args[4]

pcaobj_file <- file.path(out_dir, "pcair_pcaobj.RDS")
unrels_file <- file.path(out_dir, "pcair_unrelated_ids.txt")
pcrelate_file <- file.path(out_dir, "pcrelate_kinship.RDS")
coords_file <- file.path(out_dir, "pcair_coordinates.tsv")
plot_out <- file.path(out_dir, "figures", "pcair_pcs.svg")

suppressPackageStartupMessages({
    library(SNPRelate)
    library(gdsfmt)
    library(SeqArray)
    library(GENESIS)
    library(SeqVarTools)
    library(BiocParallel)
    library(tidyverse)
})

cat("Converting PLINK to GDS format...\n")
snpgdsBED2GDS(
    bed.fn = paste0(bed_prefix, ".bed"),
    bim.fn = paste0(bed_prefix, ".bim"),
    fam.fn = paste0(bed_prefix, ".fam"),
    out.gdsfn = gds_file
)

cat("Running PC-AiR...\n")
genofile <- snpgdsOpen(gds_file)

kingmat <- snpgdsIBDKING(genofile, verbose = TRUE)
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
kingmat2 <- kingToMatrix(kingmat, sample.id)

pcaobj <- pcair(genofile, kinobj = kingmat2, divobj = NULL, verbose = TRUE)

saveRDS(pcaobj, file = pcaobj_file)
write.table(pcaobj$unrels,
            file = unrels_file,
            quote = FALSE, row.names = FALSE, col.names = FALSE)

snpgdsClose(genofile)

cat("Converting to SeqArray format for PC-Relate...\n")
seqSNP2GDS(gds.fn = gds_file, out.fn = seq_gds_file, verbose = TRUE)

cat("Running PC-Relate...\n")
gdsfmt::showfile.gds(closeall = TRUE)

seqfile <- seqOpen(seq_gds_file, allow.duplicate = TRUE)
on.exit({ try(seqClose(seqfile), silent = TRUE) }, add = TRUE)

seqResetFilter(seqfile)

pcair_obj <- readRDS(pcaobj_file)
pcs_all <- pcair_obj$eigenvectors %||% pcair_obj$vectors

samp_active <- seqGetData(seqfile, "sample.id")
match_idx <- match(samp_active, rownames(pcs_all))
pcs_mat <- pcs_all[match_idx, , drop = FALSE]

unrels <- pcair_obj$unrels %||% attr(pcair_obj, "unrel.set")
training_set <- intersect(unrels, samp_active)

seqData <- SeqVarData(seqfile)
seqIter <- SeqVarBlockIterator(seqData, variantBlock = 20000L)

BPPARAM <- BiocParallel::SerialParam()

relate <- pcrelate(
    seqIter,
    pcs = pcs_mat,
    training.set = training_set,
    ibd.probs = FALSE,
    scale = "variant",
    small.samp.correct = TRUE,
    BPPARAM = BPPARAM,
    verbose = TRUE
)

saveRDS(relate, file = pcrelate_file)

cat("Saving PC coordinates and plotting...\n")
pcair_obj <- readRDS(pcaobj_file)
pcs <- pcair_obj$eigenvectors

n_pcs <- min(10, ncol(pcs))
pc_names <- paste0("PC", 1:n_pcs)

coords <- as_tibble(pcs[, 1:n_pcs], rownames = "IID")
colnames(coords)[-1] <- pc_names

write_tsv(coords, coords_file)

cat("PC-AiR and PC-Relate completed successfully\n")
