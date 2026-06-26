#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  cat("Usage: Rscript scripts/join_bim_rsid.R <reference.bim> <target.bim> <output.bim>\n", file = stderr())
  quit(status = 1)
}

ref_path <- args[1]
tgt_path <- args[2]
out_path <- args[3]

# Read all as character to handle both BIM (CHR ID CM POS REF ALT)
# and PVAR (#CHROM POS ID REF ALT) formats.
ref <- read.table(ref_path, header = FALSE, sep = "\t",
                  colClasses = "character", comment.char = "")
tgt <- read.table(tgt_path, header = FALSE, sep = "\t",
                  colClasses = "character", comment.char = "")
# Strip any comment/header rows (PVAR starts with #CHROM)
ref <- ref[!grepl("^#", ref[[1]]), ]
tgt <- tgt[!grepl("^#", tgt[[1]]), ]

# Detect format: BIM (CHR ID CM POS REF ALT, 6 cols) vs PVAR (#CHROM POS ID REF ALT, 5 cols)
for (d in c("ref", "tgt")) {
  dat <- get(d)
  if (ncol(dat) == 6) {
    names(dat) <- c("chr", "id", "cm", "pos", "ref", "alt")
  } else {
    names(dat) <- c("chr", "pos", "id", "ref", "alt")
    dat$cm <- "0"
  }
  assign(d, dat)
}

ref$chr_norm <- sub("^chr", "", ref$chr, ignore.case = TRUE)
tgt$chr_norm <- sub("^chr", "", tgt$chr, ignore.case = TRUE)

ref$exact_key <- paste(ref$chr_norm, ref$pos, toupper(ref$ref), toupper(ref$alt), sep = ":")
ref$pos_key  <- paste(ref$chr_norm, ref$pos, sep = ":")
tgt$exact_key <- paste(tgt$chr_norm, tgt$pos, toupper(tgt$ref), toupper(tgt$alt), sep = ":")
tgt$pos_key  <- paste(tgt$chr_norm, tgt$pos, sep = ":")

# Multi-allelic positions in ref (ambiguous for position-based matching)
ambig_pos <- unique(ref$pos_key[duplicated(ref$pos_key)])

# Lookup: exact key → RSID
exact_map <- setNames(ref$id, ref$exact_key)

# Lookup: pos key → RSID (unambiguous positions only)
unambig <- ref[!ref$pos_key %in% ambig_pos, ]
pos_map <- setNames(unambig$id, unambig$pos_key)

# 1. Exact match
tgt$rsid <- exact_map[tgt$exact_key]

# 2. Position match (fallback for unmatched, unambiguous positions only)
no_exact <- is.na(tgt$rsid)
tgt$rsid[no_exact] <- pos_map[tgt$pos_key[no_exact]]

# Replace unmatched IDs with original (keeps chr:pos:ref:alt or whatever was there)
tgt$rsid[is.na(tgt$rsid)] <- tgt$id[is.na(tgt$rsid)]

# Write output
write.table(tgt[, c("chr", "rsid", "cm", "pos", "ref", "alt")],
            out_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Report
total <- nrow(tgt)
matched_exact <- sum(!is.na(exact_map[tgt$exact_key]))
assigned_all  <- sum(tgt$rsid != tgt$id)
matched_pos   <- assigned_all - matched_exact
unmatched     <- total - assigned_all

cat(sprintf("Total target variants:  %d\n", total), file = stderr())
cat(sprintf("  Exact chr:pos:ref:alt: %d\n", matched_exact), file = stderr())
cat(sprintf("  Position only:         %d\n", matched_pos), file = stderr())
cat(sprintf("  Unmatched:             %d\n", unmatched), file = stderr())
cat(sprintf("  RSIDs assigned:        %d\n", matched_exact + matched_pos), file = stderr())
