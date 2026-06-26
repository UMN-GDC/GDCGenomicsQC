#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  cat("Usage: Rscript scripts/join_bim_rsid.R <dbsnp.txt.gz> <target.bim|pvar> <output.bim>\n", file = stderr())
  quit(status = 1)
}

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("data.table package is required. Install with install.packages('data.table')", call. = FALSE)
}

ref_path <- args[1]
tgt_path <- args[2]
out_path <- args[3]

# ═══════════════════════════════════════════════════════════════════
# 1. Read dbSNP reference (txt.gz, UCSC dbSNP table format)
# ═══════════════════════════════════════════════════════════════════
# Only keep columns 2 (chrom), 3 (chromStart, 0-based), 5 (rsID), 10 (observed)
cat("Reading dbSNP reference ... ", file = stderr())
ref <- data.table::fread(ref_path, sep = "\t", header = FALSE,
                         select = c(2, 3, 5, 10),
                         col.names = c("chrom", "pos0", "rsid", "observed"),
                         quote = "")
cat(sprintf("%d rows\n", nrow(ref)), file = stderr())

# 0-based → 1-based
ref$pos <- ref$pos0 + 1L
ref$chr_norm <- sub("^chr", "", ref$chrom, ignore.case = TRUE)

# Parse observed = "A/G" or "-/C" or "A/C/G" (skip multi-allelic)
has_slash <- grepl("/", ref$observed)
ref <- ref[has_slash, ]
allele1 <- toupper(sub("/.*", "", ref$observed))
allele2 <- toupper(sub(".*/", "", ref$observed))
# Keep only exactly 2 distinct alleles
biallelic <- allele1 != allele2 & !grepl("/", allele1) & !grepl("/", allele2)
ref <- ref[biallelic, ]
allele1 <- allele1[biallelic]
allele2 <- allele2[biallelic]

# Build lookup with both orientations (ref→alt and alt→ref)
cat("Building lookup table ... ", file = stderr())
ref_fwd <- data.frame(chr_norm = ref$chr_norm, pos = ref$pos,
                       ref = allele1, alt = allele2, rsid = ref$rsid,
                       stringsAsFactors = FALSE)
ref_rev <- data.frame(chr_norm = ref$chr_norm, pos = ref$pos,
                       ref = allele2, alt = allele1, rsid = ref$rsid,
                       stringsAsFactors = FALSE)
ref <- rbind(ref_fwd, ref_rev)
rm(ref_fwd, ref_rev, allele1, allele2, biallelic, has_slash)

ref$exact_key <- paste(ref$chr_norm, ref$pos, ref$ref, ref$alt, sep = ":")
ref$pos_key  <- paste(ref$chr_norm, ref$pos, sep = ":")
cat(sprintf("%d lookup entries\n", nrow(ref)), file = stderr())

# ═══════════════════════════════════════════════════════════════════
# 2. Read target (BIM or PVAR)
# ═══════════════════════════════════════════════════════════════════
find_header_line <- function(path) {
  lines <- readLines(path, n = 100, warn = FALSE)
  for (i in seq_along(lines)) {
    if (grepl("^#CHROM", lines[i])) return(i)
  }
  0L
}
tgt_hdr <- find_header_line(tgt_path)
is_tgt_pvar <- tgt_hdr > 0

tgt <- read.table(tgt_path, header = FALSE, sep = "\t",
                  colClasses = "character", comment.char = "",
                  skip = tgt_hdr)

if (is_tgt_pvar) {
  names(tgt) <- c("chr", "pos", "id", "ref", "alt", "cm")[1:ncol(tgt)]
} else {
  names(tgt) <- c("chr", "id", "cm", "pos", "ref", "alt")
}

tgt$chr_norm <- sub("^chr", "", tgt$chr, ignore.case = TRUE)
tgt$exact_key <- paste(tgt$chr_norm, tgt$pos, toupper(tgt$ref), toupper(tgt$alt), sep = ":")
tgt$pos_key  <- paste(tgt$chr_norm, tgt$pos, sep = ":")

# ═══════════════════════════════════════════════════════════════════
# 3. Join
# ═══════════════════════════════════════════════════════════════════
ambig_pos <- unique(ref$pos_key[duplicated(ref$pos_key)])

exact_map <- setNames(ref$rsid, ref$exact_key)
unambig   <- ref[!ref$pos_key %in% ambig_pos, ]
pos_map   <- setNames(unambig$rsid, unambig$pos_key)

tgt$rsid <- exact_map[tgt$exact_key]

no_exact <- is.na(tgt$rsid)
tgt$rsid[no_exact] <- pos_map[tgt$pos_key[no_exact]]

tgt$rsid[is.na(tgt$rsid)] <- tgt$id[is.na(tgt$rsid)]

# ═══════════════════════════════════════════════════════════════════
# 4. Write output
# ═══════════════════════════════════════════════════════════════════
write.table(tgt[, c("chr", "rsid", "cm", "pos", "ref", "alt")],
            out_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# ═══════════════════════════════════════════════════════════════════
# 5. Report
# ═══════════════════════════════════════════════════════════════════
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
