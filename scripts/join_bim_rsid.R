#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  cat("Usage: Rscript scripts/join_bim_rsid.R <reference.bim> <target.bim> <output.bim>\n", file = stderr())
  quit(status = 1)
}

ref_path <- args[1]
tgt_path <- args[2]
out_path <- args[3]

# Find the #CHROM header line number so we can skip all preceding ## lines
find_header_line <- function(path) {
  lines <- readLines(path, n = 100, warn = FALSE)
  for (i in seq_along(lines)) {
    if (grepl("^#CHROM", lines[i])) return(i)
  }
  0L  # BIM — no header
}
ref_hdr <- find_header_line(ref_path)
tgt_hdr <- find_header_line(tgt_path)

is_ref_pvar <- ref_hdr > 0
is_tgt_pvar <- tgt_hdr > 0

# For PVAR, skip past #CHROM (all preceding ## lines are ignored too).
# For BIM, skip=0 reads from line 1.
ref <- read.table(ref_path, header = FALSE, sep = "\t",
                  colClasses = "character", comment.char = "",
                  skip = ref_hdr)
tgt <- read.table(tgt_path, header = FALSE, sep = "\t",
                  colClasses = "character", comment.char = "",
                  skip = tgt_hdr)

# Assign column names based on detected format
for (d in c("ref", "tgt")) {
  is_pvar <- if (d == "ref") is_ref_pvar else is_tgt_pvar
  dat <- get(d)
  if (is_pvar) {
    names(dat) <- c("chr", "pos", "id", "ref", "alt", "cm")[1:ncol(dat)]
  } else {
    names(dat) <- c("chr", "id", "cm", "pos", "ref", "alt")
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
