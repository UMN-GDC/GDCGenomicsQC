library(tidyverse)

dir <- "/scratch.global/coffm049/r25Pipeline"
name <- "full"

work_dir <- paste0(dir, "/02-localAncestry")
setwd(work_dir)

print("Locating RFMix .lai.msp.tsv and .lai.fb.tsv files...")
msp_files <- list.files(pattern = "\\.lai\\.msp\\.tsv$")
fb_files <- list.files(pattern = "\\.lai\\.fb\\.tsv$")
if(length(msp_files) == 0) stop("No .lai.msp.tsv files found in directory!")
if(length(fb_files) == 0) stop("No .lai.fb.tsv files found in directory!")

print(paste("Found", length(msp_files), "chromosome files"))

# Get populations from fb header
fb_header <- readLines(fb_files[1], n = 1)
populations <- sub("#reference_panel_population: ", "", fb_header)
populations <- unlist(str_split(populations, "\\s+"))
populations <- populations[populations != ""]
print(paste("Populations:", paste(populations, collapse = ", ")))

# Use first chromosome only for testing
print("Reading first chromosome only for testing...")

# Read msp and calculate segment length
all_msp <- read_tsv(msp_files[1], skip = 1, show_col_types = FALSE) |>
  mutate(segment_length = epos - spos, epos = epos-1) |>
  select(chr = `#chm`, spos, epos, segment_length) |>
  distinct()

# Read fb
all_fb <- read_tsv(fb_files[1], skip = 1, show_col_types = FALSE)

# For each SNP in fb, find which msp segment it falls into (spos <= physical_position <= epos)
fb_ancestry <- all_fb |>
  left_join(
    all_msp,
    join_by(chromosome == chr, between(physical_position, spos, epos))
  ) |>
  filter(!is.na(segment_length))

sample_cols_fb <- names(all_fb)[5:ncol(all_fb)]
sample_ids <- sub(":::hap.*", "", sample_cols_fb) |> unique()
print(paste("Found", length(sample_ids), "samples"))

print("Calculating length-weighted global ancestry...")

# Pivot: wide -> (sample_id x hap) per position -> average haplotypes -> weight by segment_length
fb_by_pos <- fb_ancestry |>
  select(chromosome, physical_position, segment_length, all_of(sample_cols_fb)) |>
  pivot_longer(
    cols = -c(chromosome, physical_position, segment_length),
    names_to = c("sample_id", "hap", "pop"),
    names_pattern = "^([^:]+):::hap([12]):::(.+)$"
  ) |>
  pivot_wider(
    id_cols = c(chromosome, physical_position, segment_length, sample_id, pop),
    names_from = hap,
    values_from = value,
    names_prefix = "hap"
  ) |>
  mutate(ancestry_val = (hap1 + hap2) / 2)

# Calculate total length from unique positions per sample (before population expansion)
total_length <- fb_by_pos |>
  distinct(chromosome, physical_position, sample_id, segment_length) |>
  group_by(sample_id) |>
  summarise(total_len = sum(segment_length, na.rm = TRUE), .groups = "drop")

# Calculate weighted sum per sample/pop, then divide by total length
weighted <- fb_by_pos |>
  group_by(sample_id, pop) |>
  summarise(weighted_sum = sum(segment_length * ancestry_val, na.rm = TRUE), .groups = "drop") |>
  left_join(total_length, by = "sample_id") |>
  mutate(proportion = weighted_sum / total_len) |>
  select(sample_id, pop, proportion) |>
  pivot_wider(names_from = pop, values_from = proportion, values_fill = 0) |>
  rename(IID = sample_id)

# Check sums
print("Checking row sums...")
row_sums <- weighted |> select(-IID) |> rowSums()
print(summary(row_sums))

print("First few rows:")
print(head(weighted))
