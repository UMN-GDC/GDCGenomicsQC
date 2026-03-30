library(karyoploteR)
library(tidyverse)
library(argparse)
library(vroom)

parser <- ArgumentParser(description="Plot local ancestry estimates on karyotype")
parser$add_argument("--msp-dir", "-m", required=TRUE, help="Directory with .lai.msp.tsv files")
parser$add_argument("--sample", "-s", required=TRUE, help="Sample ID")
parser$add_argument("--chromosomes", "-c", default="1-22", help="Chromosome range (e.g., 1-22)")
parser$add_argument("--output", "-o", default="karyotype.pdf", help="Output PDF file")
parser$add_argument("--genome", "-g", default="hg38", help="Genome build (hg38 or hg19)")
args <- parser$parse_args()

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

hap0_col <- paste0(args$sample, ".0")
hap1_col <- paste0(args$sample, ".1")

# Parse chromosomes
if(grepl("-", args$chromosomes)) {
  rng <- strsplit(args$chromosomes, "-")[[1]]
  chroms <- paste0("chr", as.integer(rng[1]):as.integer(rng[2]))
} else {
  chroms <- strsplit(args$chromosomes, ",")[[1]]
  chroms <- ifelse(grepl("^chr", chroms), chroms, paste0("chr", chroms))
}

# Start PDF device BEFORE plotting
pdf(args$output, width=14, height=10)

# Initialize karyotype plot
kp <- plotKaryotype(genome=args$genome, chromosomes=chroms, plot.type=2,
                    main=paste("Local Ancestry:", args$sample))

kpAddCytobandLines(kp, layerMargin=0)

# Prepare for plotting
# Panel 1 (above ideogram) is Hap 0
# Panel 2 (below ideogram) is Hap 1
y0 <- 0
y1 <- 1

# Optional: background for data panels
kpDataBackground(kp, data.panel=1, r0=y0, r1=y1, col="#F8F8F8", border=NA)
kpDataBackground(kp, data.panel=2, r0=y0, r1=y1, col="#F8F8F8", border=NA)

# Try to extract population names from MSP header if available
pop_map <- NULL
for(chr in chroms) {
  msp_file <- file.path(args$msp_dir, paste0(chr, ".lai.msp.tsv"))
  if(file.exists(msp_file)) {
    header_line <- readLines(msp_file, n=1)
    if(grepl("populations:", header_line)) {
      pops <- sub(".*populations: ", "", header_line) |> strsplit(" ") |> unlist()
      # pops looks like c("EUR=0", "AFR=1")
      pop_map <- setNames(sub("=.*", "", pops), sub(".*=", "", pops))
      break
    }
  }
}

for(chr in chroms) {
  msp_file <- file.path(args$msp_dir, paste0(chr, ".lai.msp.tsv"))
  if(!file.exists(msp_file)) {
    message("Warning: Missing MSP file for ", chr)
    next
  }
  
  # Read MSP file (skip comment line)
  # col_select uses strings for robustness
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
  # If we found pop names, use them: "EUR (0)"
  match_idx <- match(legend_labels, names(pop_map))
  legend_labels[!is.na(match_idx)] <- paste0(pop_map[match_idx[!is.na(match_idx)]], " (", legend_labels[!is.na(match_idx)], ")")
}

# Only show legend for used colors
legend("bottomright", legend=legend_labels[1:length(ancestry_colors)], 
       fill=ancestry_colors, title="Ancestry", bty="n", cex=0.8)

# Add Haplotype labels
kpAddText(kp, chr=chroms[1], x=0, y=0.5, labels="Hap 0", data.panel=1, pos=2, cex=0.8)
kpAddText(kp, chr=chroms[1], x=0, y=0.5, labels="Hap 1", data.panel=2, pos=2, cex=0.8)

dev.off()

message("Success: Plot saved to ", args$output)
