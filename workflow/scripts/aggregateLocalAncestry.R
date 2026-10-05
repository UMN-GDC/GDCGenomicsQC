library(data.table)

args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]
name <- "full"

work_dir <- paste0(dir, "/02-localAncestry")
setwd(work_dir)

print("Locating RFMix .lai.msp.tsv and .lai.fb.tsv files...")
msp_files <- list.files(pattern = "\\.lai\\.msp\\.tsv$")
fb_files <- list.files(pattern = "\\.lai\\.fb\\.tsv$")
if (length(msp_files) == 0) stop("No .lai.msp.tsv files found in directory!")
if (length(fb_files) == 0) stop("No .lai.fb.tsv files found in directory!")

print(paste("Found", length(msp_files), "chromosome files"))

# Populations from first fb header line
fb_header <- readLines(fb_files[1], n = 1)
populations <- sub("#reference_panel_population: ", "", fb_header)
populations <- unlist(strsplit(populations, "\\s+"))
populations <- populations[populations != ""]
populations <- populations[!grepl("^#", populations)]
print(paste("Populations:", paste(populations, collapse = ", ")))

# Distinct msp segments (per chromosome) with lengths
print("Reading msp files for segment boundaries...")
msp_list <- lapply(msp_files, function(f) {
  d <- fread(f, skip = 1, showProgress = FALSE)
  setnames(d, 1, "chr")
  d[, segment_length := epos - spos]
  unique(d[, .(chr, spos, epos, segment_length)])
})
msp_seg <- rbindlist(msp_list)
setkey(msp_seg, chr, spos, epos)

accum <- list()  # sample -> named pop sums
total_len <- list()

for (f in fb_files) {
  print(paste("Processing", f, "..."))
  fb <- fread(f, skip = 1, showProgress = FALSE)
  fb_cols <- names(fb)
  # Map each fb row to its msp segment length (non-equi join)
  fb[, seglen := msp_seg[.SD, on = .(chr = chromosome, spos <= physical_position, epos >= physical_position),
                         x.segment_length]]
  fb <- fb[!is.na(seglen)]
  if (nrow(fb) == 0) next
  w <- fb$seglen
  tot_block <- sum(w)
  file_sums <- list()  # pop -> named sample vector for this file
  for (pop in populations) {
    h1_cols <- grep(paste0(":::hap1:::", pop, "$"), fb_cols, value = TRUE)
    h2_cols <- grep(paste0(":::hap2:::", pop, "$"), fb_cols, value = TRUE)
    if (length(h1_cols) == 0) next
    m1 <- as.matrix(fb[, ..h1_cols])
    m2 <- as.matrix(fb[, ..h2_cols])
    storage.mode(m1) <- "double"
    storage.mode(m2) <- "double"
    a <- (m1 + m2) / 2
    a[is.na(a)] <- 0
    ws <- colSums(a * w)
    names(ws) <- sub(":::hap.*", "", names(ws))
    file_sums[[pop]] <- ws
    rm(m1, m2, a, ws)
    gc()
  }
  # Merge this file's sums into global accumulators (totals added once)
  all_samp <- unique(unlist(lapply(file_sums, names)))
  for (s in all_samp) {
    if (is.null(accum[[s]])) accum[[s]] <- setNames(rep(0, length(populations)), populations)
    for (pop in names(file_sums)) {
      v <- file_sums[[pop]][s]
      if (!is.na(v)) accum[[s]][pop] <- accum[[s]][pop] + v
    }
    total_len[[s]] <- (if (is.null(total_len[[s]])) 0 else total_len[[s]]) + tot_block
  }
  rm(fb)
  gc()
}

sample_ids <- names(accum)
print(paste("Found", length(sample_ids), "samples"))

print("Calculating length-weighted global ancestry...")
mat <- t(vapply(sample_ids, function(s) accum[[s]] / total_len[[s]], numeric(length(populations))))
colnames(mat) <- populations
ancestry_df <- data.frame(IID = sample_ids, mat, check.names = FALSE)

output_mat_path <- paste0(dir, "/02-localAncestry/ancestry_", name, ".txt")
write.table(ancestry_df, file = output_mat_path, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

print("RFMix Global Ancestry calculation complete!")
