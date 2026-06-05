# Consolidation of gai.R / gai2.R / gai3.R
#
# Converts RFMix output (.msp.tsv or .rfmix.Q) to per-sample ancestry estimates
# and writes study.{name}.unrelated.comm.popu
#
# Usage: Rscript rfmix_to_ancestry.R <dir> <name> [method]
#   method: "segment" (default, from gai.R) or "qmean" (from gai2.R/gai3.R)
#   "segment" uses segment-length-weighted ancestry proportions from .msp.tsv
#   "qmean"   uses mean ancestry across chromosomes from .rfmix.Q

library(dplyr)
library(magrittr)
library(tidyr)
library(purrr)
library(readr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]
name <- args[2]
method <- ifelse(length(args) >= 3, args[3], "segment")

work_dir <- paste0(dir, "/rfmix")
setwd(work_dir)

if (method == "segment") {
  # ---- gai.R style: segment-length-weighted from .msp.tsv ----
  msp_list <- list.files(pattern = ".msp.tsv")
  ancestry_code <- system(paste("head -n 1 ", msp_list[1]), intern = TRUE)
  ancestry_code <- sub(".*:", "", ancestry_code)
  ancestry_code <- gsub("\t", ",", ancestry_code)
  ancestry_code <- unlist(strsplit(ancestry_code, ","))
  ancestry_code <- gsub("\\s+", "", ancestry_code)
  code <- sub("=.*", "", ancestry_code)
  number <- sub(".*=", "", ancestry_code)

  msp_colnames <- system(paste("head -n 2 ", msp_list[1]), intern = TRUE)
  msp_colnames <- msp_colnames[2]
  msp_colnames <- sub("#", "", msp_colnames)
  msp_colnames <- gsub("\t", ",", msp_colnames)
  msp_colnames <- unlist(strsplit(msp_colnames, ","))

  msp_data <- list.files(pattern = ".msp.tsv") %>%
    lapply(read.table, header = FALSE, sep = "\t") %>%
    bind_rows()

  colnames(msp_data) <- msp_colnames

  segment_data <- msp_data %>%
    transmute(chm, spos, epos, slength = epos - spos, sgpos, egpos, n_snps = `n snps`)
  msp_data <- msp_data %>%
    select(-c(chm, spos, epos, sl = epos - spos, sgpos, egpos, `n snps`))

  nsamples <- ncol(msp_data) / 2
  hap_mat <- matrix(0, nrow = ncol(msp_data), ncol = length(code))
  rownames(hap_mat) <- colnames(msp_data)
  colnames(hap_mat) <- code

  length_vec <- segment_data$slength
  length_tot <- sum(length_vec)
  prop_mat <- matrix(length_vec / length_tot, nrow = 1)

  for (i in seq_along(code)) {
    ancestry_array <- matrix(as.numeric(msp_data == number[i]), nrow = nrow(msp_data))
    hap_mat[, i] <- prop_mat %*% ancestry_array
  }

  ancestry_mat <- t(sapply(seq_len(nsamples), function(x) {
    (hap_mat[2 * x - 1, ] + hap_mat[2 * x, ]) / 2
  }))
  rownames(ancestry_mat) <- str_extract(rownames(hap_mat), "^[^.]+")[2 * seq_len(nsamples)]

} else {
  # ---- qmean style: mean ancestry across chromosomes from .rfmix.Q ----
  Qfiles <- list.files(pattern = ".rfmix.Q")
  qnames <- system(paste("head -n 2 ", Qfiles[1]), intern = TRUE)
  qnames <- qnames[2]
  qnames <- sub("#", "", qnames)
  qnames <- gsub("\t", ",", qnames)
  qnames <- unlist(strsplit(qnames, ","))
  n_anc <- length(qnames) - 1

  Q_list <- lapply(Qfiles, function(f) {
    dat <- read.table(f, header = FALSE)
    colnames(dat) <- qnames
    return(dat)
  })

  sample_vec <- Q_list[[1]][, 1]

  Q_array <- array(0,
    dim = c(length(sample_vec), n_anc, length(Q_list)))
  for (i in seq_along(Q_list)) {
    Q_array[, , i] <- apply(Q_list[[i]][, -1, drop = FALSE], 2, as.numeric)
  }

  Q_data <- apply(Q_array, c(1, 2), mean)
  ancestry_mat <- Q_data
  rownames(ancestry_mat) <- sample_vec
}

colnames(ancestry_mat) <- code

# Write ancestry proportions
output_path <- paste0(dir, "/ancestry_", name, ".txt")
write.table(ancestry_mat, file = output_path, row.names = TRUE, col.names = TRUE, quote = FALSE)

# Build ancestry decision table
index <- apply(ancestry_mat, 1, which.max)
ancestry_decision <- data.frame(
  ID = names(index),
  code_number = unname(index),
  ancestry = code[index],
  prediction_percentage = apply(ancestry_mat, 1, max) * 100
)
ancestry_decision$ancestry[ancestry_decision$prediction_percentage < 80] <- "Other"

fam_name <- paste0("study.", name, ".unrelated.fam")
fam_path <- paste0(dir, "/relatedness/", fam_name)
fam_file <- read.table(fam_path, header = FALSE)
colnames(fam_file) <- c("FID", "IID", "MID", "PID", "gender", "phenotype")
fam_file$ID <- paste0(fam_file$FID, "_", fam_file$IID)

joined_file <- fam_file %>%
  inner_join(ancestry_decision, by = "ID") %>%
  select(FID, IID, ancestry, prediction_percentage, gender, phenotype)

output_path <- paste0(dir, "/ancestry_estimation/study.", name, ".unrelated.comm.popu")
dir.create(dirname(output_path), showWarnings = FALSE)
write.table(joined_file, file = output_path, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
