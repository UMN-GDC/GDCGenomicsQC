#!/usr/bin/env Rscript
library(karyoploteR)
library(tidyverse)
library(argparse)

parser <- ArgumentParser(description="Plot local ancestry estimates on karyotype")
parser$add_argument("--msp-dir", "-m", required=TRUE, help="Directory with .lai.msp.tsv files")
parser$add_argument("--sample", "-s", required=TRUE, help="Sample ID")
parser$add_argument("--chromosomes", "-c", default="1-22", help="Chromosome range (e.g., 1-22)")
parser$add_argument("--output", "-o", default="karyotype.pdf", help="Output PDF file")
parser$add_argument("--genome", "-g", default="hg38", help="Genome build (hg38 or hg19)")
args <- parser$parse_args()

ancestry_colors <- c(
  "0" = "#E69F00",
  "1" = "#56B4E9",
  "2" = "#009E73",
  "3" = "#CC79A7"
)

hap0_col <- paste0(args$sample, ".0")
hap1_col <- paste0(args$sample, ".1")

chroms <- if(grepl("-", args$chromosomes)) {
  rng <- strsplit(args$chromosomes, "-")[[1]]
  paste0("chr", as.integer(rng[1]):as.integer(rng[2]))
} else {
  strsplit(args$chromosomes, ",")[[1]] |>
    (\(x) if(grepl("^chr", x[1])) x else paste0("chr", x))()
}

chr_str <- paste0("chr", 1:22)

kp <- plotKaryotype(genome=args$genome, chromosomes=chr_str, plot.type=2,
                    main=paste("Local Ancestry:", args$sample))

kpAddCytobandLines(kp, layerMargin=0)

kpDataBackground(kp, data.panel=1, r0=0.8, r1=1)
kpDataBackground(kp, data.panel=2, r0=0, r1=0.2)

for(chr in chroms) {
  msp_file <- file.path(args$msp_dir, paste0(chr, ".lai.msp.tsv"))
  if(!file.exists(msp_file)) next
  
  msp <- vroom(msp_file, show_col_types=FALSE, skip=1,
               col_select=c("#chm", spos, epos, !!hap0_col, !!hap1_col)) |>
    rename(chr="#chm", start=spos, end=epos, hap0=!!hap0_col, hap1=!!hap1_col)
  
  for(i in 1:nrow(msp)) {
    anc0 <- as.character(msp$hap0[i])
    anc1 <- as.character(msp$hap1[i])
    
    if(!is.na(anc0) && anc0 %in% names(ancestry_colors)) {
      kpRect(kp, chr=msp$chr[i], x0=msp$start[i], x1=msp$end[i],
             y0=0.8, y1=1, col=ancestry_colors[anc0], data.panel=1)
    }
    if(!is.na(anc1) && anc1 %in% names(ancestry_colors)) {
      kpRect(kp, chr=msp$chr[i], x0=msp$start[i], x1=msp$end[i],
             y0=0, y1=0.2, col=ancestry_colors[anc1], data.panel=2)
    }
  }
}

for(i in seq_along(ancestry_colors)) {
  kpRect(kp, chr="chr1", x0=0, x1=1e6, y0=0.8+i*0.03, y1=0.83+i*0.03,
         col=ancestry_colors[i], data.panel=1)
  kpAddText(kp, labels=names(ancestry_colors)[i], chr="chr1", 
            x=5e6, y=0.8+i*0.03+0.015, data.panel=1, cex=0.7)
}

pdf(args$output, width=14, height=10)
print(kp)
dev.off()

message("Saved: ", args$output)
