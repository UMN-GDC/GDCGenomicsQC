#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 7) {
    stop("Usage: Rscript run_pcair_pcrelate.R <bed_prefix> <out_dir> <color_col> <pheno_file> <plot_out> <gds_file> <seq_gds_file>\n")
}

bed_prefix <- args[1]
out_dir <- args[2]
color_col <- args[3]
pheno_file <- args[4]
plot_out <- args[5]
gds_file <- args[6]
seq_gds_file <- args[7]

pcaobj_file <- file.path(out_dir, "pcair_pcaobj.RDS")
unrels_file <- file.path(out_dir, "pcair_unrelated_ids.txt")
pcrelate_file <- file.path(out_dir, "pcrelate_kinship.RDS")
coords_file <- file.path(out_dir, "pcair_coordinates.tsv")

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

color_col <- ifelse(color_col == "None" || color_col == "", NA, color_col)
pheno_file <- ifelse(pheno_file == "None" || pheno_file == "", NA, pheno_file)

if (!is.na(color_col) && !is.na(pheno_file) && file.exists(pheno_file)) {
    pheno <- read_tsv(pheno_file)
    if (color_col %in% colnames(pheno)) {
        coords <- coords |> left_join(pheno |> select(IID, all_of(color_col)), by = "IID")
        color_var <- color_col
    } else {
        coords$color_group <- "all"
        color_var <- "color_group"
    }
} else {
    coords$color_group <- "all"
    color_var <- "color_group"
}

p <- ggplot(coords, aes(x = PC1, y = PC2, color = .data[[color_var]])) +
    geom_point(alpha = 0.6, size = 2) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(title = "PC-AiR: Internal Sample PCs")

ggsave(plot_out, plot = p, dpi = 300, width = 8, height = 6)

cat("PC-AiR and PC-Relate completed successfully\n")
