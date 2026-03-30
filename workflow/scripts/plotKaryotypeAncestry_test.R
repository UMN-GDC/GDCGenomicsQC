library(karyoploteR)
library(tidyverse)
library(vroom)

# Hardcoded paths for testing
msp_dir <- "/scratch.global/coffm049/toyPipeline/02-localAncestry"
sample_id <- "SAS57_SAS57"
output_file <- "/scratch.global/coffm049/GDCGenomicsQC/test_karyotype.pdf"
genome <- "hg38"
chromosomes <- "20-21"

# Colorblind-friendly palette
ancestry_colors <- c(
  "0" = "#E69F00", # Orange
  "1" = "#56B4E9", # Sky Blue
  "2" = "#009E73", # Bluish Green
  "3" = "#CC79A7", # Reddish Purple
  "4" = "#F0E442", # Yellow
  "5" = "#0072B2", # Blue
  "6" = "#D55E00"  # Vermillion
)

hap0_col <- paste0(sample_id, ".0")
hap1_col <- paste0(sample_id, ".1")

# Parse chromosomes
if(grepl("-", chromosomes)) {
  rng <- strsplit(chromosomes, "-")[[1]]
  chroms <- paste0("chr", as.integer(rng[1]):as.integer(rng[2]))
} else {
  chroms <- strsplit(chromosomes, ",")[[1]]
  chroms <- ifelse(grepl("^chr", chroms), chroms, paste0("chr", chroms))
}

# Start PDF device BEFORE plotting
pdf(output_file, width=14, height=10)

# Initialize karyotype plot
kp <- plotKaryotype(genome=genome, chromosomes=chroms, plot.type=2,
                    main=paste("Local Ancestry:", sample_id))

# kpAddCytobandLines(kp, layerMargin=0)

# Prepare for plotting
y0 <- 0
y1 <- 1

# Optional: background for data panels
kpDataBackground(kp, data.panel=1, r0=y0, r1=y1, col="#F8F8F8", border=NA)
kpDataBackground(kp, data.panel=2, r0=y0, r1=y1, col="#F8F8F8", border=NA)

# Try to extract population names from MSP header if available
pop_map <- NULL
for(chr in chroms) {
  msp_file <- file.path(msp_dir, paste0(chr, ".lai.msp.tsv"))
  if(file.exists(msp_file)) {
    header_line <- readLines(msp_file, n=1)
    if(grepl("populations:", header_line)) {
      pops <- sub(".*populations: ", "", header_line) |> strsplit(" ") |> unlist()
      pop_map <- setNames(sub("=.*", "", pops), sub(".*=", "", pops))
      break
    }
  }
}

for(chr in chroms) {
  msp_file <- file.path(msp_dir, paste0(chr, ".lai.msp.tsv"))
  if(!file.exists(msp_file)) {
    message("Warning: Missing MSP file for ", chr)
    next
  }
  
  # Read MSP file (skip comment line)
  # msp <- head(names(vroom(msp_file, show_col_types=FALSE, skip=1, n_max = 0)), n= 30)
  msp <- vroom(msp_file, show_col_types=FALSE, skip=1,
               col_select=c("#chm", "spos", "epos", all_of(c(hap0_col, hap1_col)))) %>%
    rename(chr="#chm", start="spos", end="epos", hap0=all_of(hap0_col), hap1=all_of(hap1_col)) %>%
    mutate(chr = ifelse(grepl("^chr", chr), chr, paste0("chr", chr))) %>%
    filter(chr == !!chr)
  
  if(nrow(msp) == 0) next
  
  # Vectorized plotting per ancestry
  present_ancs <- unique(c(msp$hap0, msp$hap1))
  for(anc in present_ancs) {
    anc_str <- as.character(anc)
    if(!(anc_str %in% names(ancestry_colors))) next
    
    col <- ancestry_colors[anc_str]
    
    # Hap 0 (Data Panel 1)
    d0 <- msp[msp$hap0 == anc, ]
    if(nrow(d0) > 0) {
      kpRect(kp, chr=d0$chr, x0=d0$start, x1=d0$end, y0=y0, y1=y1, 
             col=col, border=NA, data.panel=1)
    }
    
    # Hap 1 (Data Panel 2)
    d1 <- msp[msp$hap1 == anc, ]
    if(nrow(d1) > 0) {
      kpRect(kp, chr=d1$chr, x0=d1$start, x1=d1$end, y0=y0, y1=y1, 
             col=col, border=NA, data.panel=2)
    }
  }
}

# Add legend
legend_labels <- names(ancestry_colors)
if(!is.null(pop_map)) {
  match_idx <- match(legend_labels, names(pop_map))
  legend_labels[!is.na(match_idx)] <- paste0(pop_map[match_idx[!is.na(match_idx)]], " (", legend_labels[!is.na(match_idx)], ")")
}

legend("bottomright", legend=legend_labels[1:length(ancestry_colors)], 
       fill=ancestry_colors, title="Ancestry", bty="n", cex=0.8)

# Add Haplotype labels
# kpAddText(kp, chr=chroms[1], x=0, y=0.5, labels="Hap 0", data.panel=1, pos=2, cex=0.8)
# kpAddText(kp, chr=chroms[1], x=0, y=0.5, labels="Hap 1", data.panel=2, pos=2, cex=0.8)

dev.off()

message("Success: Test plot saved to ", output_file)
