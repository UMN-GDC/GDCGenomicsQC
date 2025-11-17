#!/usr/bin/env Rscript

## PC-Relate pipeline with LD pruning (SNPRelate + SeqArray + GENESIS)
##
## Usage (from your SLURM script):
##   Rscript pca_relate_pipeline.R ${WORK} ${NAME}
##
## Expects (from pca_ir_pipeline.R):
##   WORK/gds/NAME.gds           : SNPRelate GDS
##   WORK/gds/NAME_seq.gds       : SeqArray GDS
##   WORK/pca_ir/NAME_pcaobj.RDS : PC-AiR object
##
## Output:
##   WORK/pca_ir/NAME_pcrelate.RDS : PC-Relate result

args <- commandArgs(trailingOnly = TRUE)
WORK <- args[1]
NAME <- args[2]

suppressPackageStartupMessages({
  library(gdsfmt)
  library(SeqArray)
  library(SeqVarTools)
  library(GENESIS)
  library(BiocParallel)
})

`%||%` <- function(x, y) if (!is.null(x)) x else y

## ============================================================
## Step 0: Clean up any open GDS handles
## ============================================================
gdsfmt::showfile.gds(closeall = TRUE)

## ---------------- Paths & parameters ----------------
seq_gds_file <- file.path(WORK, "gds", paste0(NAME, "_seq.gds"))
out_dir      <- file.path(WORK, "pca_ir")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pcair_rds    <- file.path(out_dir, paste0(NAME, "_pcaobj.RDS"))
if (!file.exists(seq_gds_file)) {
  stop("SeqArray GDS not found: ", seq_gds_file)
}
if (!file.exists(pcair_rds)) {
  stop("PC-AiR object not found: ", pcair_rds)
}

## PC-Relate knobs
variant_block <- 20000L     # variants per block in iterator

## Optional region filter (MHC on GRCh37-ish)
drop_MHC  <- TRUE
mhc_chr   <- 6
mhc_start <- 25e6
mhc_end   <- 34e6

## Helper to normalize chromosome labels
norm_chr <- function(x) {
  x <- as.character(x)
  x[x %in% c("23","X")]  <- "X"
  x[x %in% c("24","Y")]  <- "Y"
  x[x %in% c("25","MT","M","MTDNA","MTdna","Mt")] <- "MT"
  x
}

## ============================================================
## Step 1: Open SeqArray GDS and basic diagnostics
## ============================================================
cat("Opening SeqArray GDS for PC-Relate: ", seq_gds_file, "\n", sep = "")
seqfile <- SeqArray::seqOpen(seq_gds_file, allow.duplicate = TRUE)
on.exit({ try(SeqArray::seqClose(seqfile), silent = TRUE) }, add = TRUE)

## Start from a clean slate in this handle: ALL samples, ALL variants
SeqArray::seqResetFilter(seqfile)

all_samples <- SeqArray::seqGetData(seqfile, "sample.id")
all_chr     <- SeqArray::seqGetData(seqfile, "chromosome")
all_pos     <- SeqArray::seqGetData(seqfile, "position")

cat("SeqArray contains ", length(all_samples), " samples and ",
    length(all_chr), " variants.\n", sep = "")
cat("Example sample.id[1:5]: ",
    paste(head(all_samples, 5), collapse = ", "), "\n", sep = "")
cat("Example chr[1:5]: ",
    paste(head(all_chr, 5), collapse = ", "), "\n", sep = "")

if (length(all_samples) == 0L) stop("GDS has 0 samples in sample.id.")
if (length(all_chr)     == 0L) stop("GDS has 0 variants in chromosome node.")

## ============================================================
## Step 2: Variant selection via logical mask (MHC drop only)
## ============================================================
keep_mask <- rep(TRUE, length(all_chr))

if (drop_MHC) {
  cat("Dropping MHC region: chr", mhc_chr, ":", mhc_start, "-", mhc_end, "\n", sep = "")
  chr_norm <- norm_chr(all_chr)
  mhc_mask <- (chr_norm == norm_chr(mhc_chr)) &
    all_pos >= mhc_start & all_pos <= mhc_end
  keep_mask <- !mhc_mask
  
  cat("Variants remaining after MHC drop: ",
      sum(keep_mask), " / ", length(all_chr), "\n", sep = "")
}

if (sum(keep_mask) == 0L) {
  stop("Variant selection (e.g., MHC drop) left 0 variants. Loosen filters.")
}

cat("Number of variants selected (post MHC logic): ",
    sum(keep_mask), "\n", sep = "")

## ============================================================
## Step 3: Apply filters in GDS using logical variant.sel
## ============================================================
SeqArray::seqResetFilter(seqfile)  # back to ALL samples, ALL variants

cat("Setting GDS filter to all samples and ",
    sum(keep_mask), " variants via logical mask (no LD pruning)...\n", sep = "")

SeqArray::seqSetFilter(
  seqfile,
  variant.sel = keep_mask,
  verbose     = FALSE
)

## Now count active samples/variants via data
samp_active <- SeqArray::seqGetData(seqfile, "sample.id")
chr_active  <- SeqArray::seqGetData(seqfile, "chromosome")

nsamp_active <- length(samp_active)
nvar_active  <- length(chr_active)

cat("After applying variant mask: ",
    nsamp_active, " samples, ",
    nvar_active, " variants.\n", sep = "")

if (nsamp_active == 0L || nvar_active == 0L) {
  cat("DEBUG: sum(keep_mask)        = ", sum(keep_mask), "\n", sep = "")
  cat("DEBUG: length(all_samples)   = ", length(all_samples), "\n", sep = "")
  cat("DEBUG: length(all_chr)       = ", length(all_chr), "\n", sep = "")
  stop("After applying logical variant.sel filter, 0 samples or 0 variants are active. ",
       "This indicates a mismatch between keep_mask and chromosome length.")
}

## ============================================================
## Step 4: Load PC-AiR object and align PCs to SeqArray samples
## ============================================================
cat("Loading PC-AiR object from: ", pcair_rds, "\n", sep = "")
pcair_obj <- readRDS(pcair_rds)

## Try to find eigenvectors matrix in the PC-AiR object
pcs_all <- pcair_obj$eigenvectors %||% pcair_obj$vectors
if (is.null(pcs_all)) {
  stop("Could not find eigenvectors in PC-AiR object. Tried $eigenvectors and $vectors.")
}

if (is.null(rownames(pcs_all))) {
  stop("PC-AiR eigenvector matrix has no rownames; cannot align to sample.id.")
}

## Unrelated set: usually pcair_obj$unrels or attr(, 'unrel.set')
unrels <- pcair_obj$unrels %||% attr(pcair_obj, "unrel.set")
if (is.null(unrels)) {
  warning("No unrelated set found in PC-AiR object; using all samples as training.set.")
  unrels <- rownames(pcs_all)
}

## Align PCs to the current sample order in the GDS
cat("Aligning PCs to GDS sample order...\n")
match_idx <- match(samp_active, rownames(pcs_all))

if (any(is.na(match_idx))) {
  missing_ids <- samp_active[is.na(match_idx)]
  cat("First few unmatched sample IDs in PC-AiR object:\n")
  print(head(missing_ids))
  stop("Some GDS sample.id values do not appear in the PC-AiR PCs rownames; cannot align PCs.")
}

pcs_mat <- pcs_all[match_idx, , drop = FALSE]

## Training set = unrelateds that are also in the active samples
training_set <- intersect(unrels, samp_active)
if (length(training_set) == 0L) {
  warning("No overlap between PC-AiR unrelated set and active GDS samples; ",
          "using all samples as training.set.")
  training_set <- samp_active
}

cat("PC matrix has ", nrow(pcs_mat), " rows (samples) and ",
    ncol(pcs_mat), " columns (PCs).\n", sep = "")
cat("Training set size (unrelateds ∩ active samples): ",
    length(training_set), "\n", sep = "")

## ============================================================
## Step 5: Build SeqVarData + block iterator
## ============================================================
seqData <- SeqVarTools::SeqVarData(seqfile)
seqIter <- SeqVarTools::SeqVarBlockIterator(seqData, variantBlock = variant_block)

cat("Iterator will traverse ", nvar_active, " variants in blocks of ",
    variant_block, " across ", nsamp_active, " samples.\n", sep = "")
if (nvar_active <= 1000) {
  stop("Iterator sees only ", nvar_active,
       " variants — extremely low for PC-Relate.")
}

## ============================================================
## Step 6: PC-Relate (with PCs from PC-AiR, SERIAL)
## ============================================================
## IMPORTANT: use SerialParam to avoid GDS handle issues in forked workers
BPPARAM <- BiocParallel::SerialParam()

cat("Running PC-Relate on ", nsamp_active, " samples and ",
    nvar_active, " variants using 1 CPU core (SerialParam, no LD prune, MHC ",
    ifelse(drop_MHC, "dropped", "kept"), ")...\n", sep = "")

relate <- GENESIS::pcrelate(
  seqIter,
  pcs               = pcs_mat,
  training.set      = training_set,
  ibd.probs         = FALSE,
  scale             = "variant",
  small.samp.correct= TRUE,
  BPPARAM           = BPPARAM,
  verbose           = TRUE
)

## ============================================================
## Step 7: Save and close
## ============================================================
out_rds <- file.path(out_dir, paste0(NAME, "_pcrelate_seq_withPC_noLD_serial.RDS"))
saveRDS(relate, file = out_rds)
cat("Saved PC-Relate results (SeqArray + PC-AiR PCs, no LD prune, MHC ",
    ifelse(drop_MHC, "dropped", "kept"), ", SerialParam) to: ",
    out_rds, "\n", sep = "")

cat("PC-Relate pipeline finished successfully.\n")
