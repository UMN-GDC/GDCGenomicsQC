#!/usr/bin/env Rscript
library(karyoploteR)
library(tidyverse)
library(vroom)

# ============================================================================
# HARDCODE PATHS HERE
# ============================================================================
MSP_DIR <- "/scratch.global/coffm049/toyPipeline/02-localAncestry"
SAMPLE  <- "AFR57_AFR57"  # Change to your sample ID (format: ID_ID)
OUTPUT  <- "/scratch.global/coffm049/toyPipeline/04-plots/test_karyotype.pdf"
GENOME  <- "hg38"      # or "hg19"
CHROMOSOMES <- 20:21   # Use 1:22 for all chromosomes
# ============================================================================

dir.create(dirname(OUTPUT), showWarnings=FALSE, recursive=TRUE)

ancestry_colors <- c(
  "0" = "#E69F00",
  "1" = "#56B4E9",
  "2" = "#009E73",
  "3" = "#CC79A7"
)

ancestry_labels <- c(
  "0" = "EUR",
  "1" = "AFR",
  "2" = "NAT",
  "3" = "SAS"
)

hap0_col <- paste0(SAMPLE, ".0")
hap1_col <- paste0(SAMPLE, ".1")
chroms <- paste0("chr", CHROMOSOMES)
chr_str <- paste0("chr", 1:22)

message("Loading karyotype...")
kp <- plotKaryotype(genome=GENOME, chromosomes=chr_str, plot.type=2,
                    main=paste("Local Ancestry:", SAMPLE))

kpAddCytobandLines(kp, layerMargin=0)
kpDataBackground(kp, data.panel=1, r0=0.8, r1=1)
kpDataBackground(kp, data.panel=2, r0=0, r1=0.2)

message("Processing chromosomes...")
for(chr in chroms) {
  msp_file <- file.path(MSP_DIR, paste0(chr, ".lai.msp.tsv"))
  if(!file.exists(msp_file)) {
    warning("File not found: ", msp_file)
    next
  }
  
  message("  ", chr)
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

message("Adding legend...")
y_pos <- 0.82
for(i in seq_along(ancestry_colors)) {
  kpRect(kp, chr="chr1", x0=0, x1=5e6, y0=y_pos, y1=y_pos+0.025,
         col=ancestry_colors[i], data.panel=1)
  kpAddText(kp, labels=ancestry_labels[names(ancestry_colors)[i]], chr="chr1",
            x=8e6, y=y_pos+0.012, data.panel=1, cex=0.6, color="black")
  y_pos <- y_pos + 0.035
}

kpAddLabels(kp, labels="Hap1", side="left", data.panel=1, r0=0.8, r1=1, cex=0.6)
kpAddLabels(kp, labels="Hap2", side="left", data.panel=2, r0=0, r1=0.2, cex=0.6)

message("Saving to ", OUTPUT, "...")
pdf(OUTPUT, width=14, height=10)
print(kp)
dev.off()

message("Done! Saved: ", OUTPUT)
