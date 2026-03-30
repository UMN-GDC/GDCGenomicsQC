#!/usr/bin/env Rscript
library(GenomicRanges)
library(tidyverse)
library(vroom)
library(RIdeogram)
library(patchwork)
library(argparse)

parser <- ArgumentParser(description="Plot local ancestry estimates from .lai.msp.tsv files")
parser$add_argument("--input-dir", "-i", required=TRUE, help="Directory containing .lai.msp.tsv files")
parser$add_argument("--output", "-o", default="local_ancestry_plot.png", help="Output file path")
parser$add_argument("--subject", "-s", required=TRUE, help="Subject ID (e.g., SAS34825_SAS34825)")
parser$add_argument("--chromosomes", "-c", default="1-22", help="Chromosome range (e.g., 1-22 or chr1,chr2)")
parser$add_argument("--genome", "-g", default="hg38", help="Genome build (hg38 or hg19)")
args <- parser$parse_args()

genomeBuild <- args$genome
mspDir <- args$input_dir
subj <- args$subject
hap0 <- paste0(subj, ".0")
hap1 <- paste0(subj, ".1")

data(human_karyotype, package="RIdeogram")

colors <- c(
  "0" = "#E69F00",
  "1" = "#56B4E9",
  "2" = "#009E73",
  "3" = "#009E73"
)

chroms <- if(grepl("-", args$chromosomes)) {
  rng <- strsplit(args$chromosomes, "-")[[1]]
  as.integer(rng[1]):as.integer(rng[2])
} else {
  strsplit(args$chromosomes, ",")[[1]]
}

plots <- list()

for(chr in chroms) {
  chr_label <- paste0("chr", chr)
  mspPath <- file.path(mspDir, paste0(chr_label, ".lai.msp.tsv"))
  
  if(!file.exists(mspPath)) {
    warning(paste("File not found:", mspPath))
    next
  }
  
  msp <- vroom(mspPath, show_col_types = FALSE,
    col_names = TRUE, skip = 1,
    col_sel = c(`#chm`, spos, epos, sgpos, egpos,
    !!hap0, !!hap1)) |>
    rename(Chr = `#chm`, Start = spos, End = epos, Value = !!hap0, hap2 = !!hap1)
  
  p <- ideogram(karyotype = human_karyotype, overlaid = msp, 
                label_type = "marker", colorset1 = colors)
  
  svg_file <- tempfile(fileext = ".svg")
  svg_code <- capture.output(cat(p))
  writeLines(svg_code, svg_file)
  
  gg <- png::readPNG(svg_file)
  plots[[chr]] <- cowplot::ggdraw() + cowplot::background_image(gg) + cowplot::draw_plot_label(chr_label, size=16)
  
  unlink(svg_file)
}

if(length(plots) == 0) {
  stop("No plots generated. Check input files exist.")
}

combined <- wrap_plots(plots, ncol = 4) + 
  plot_annotation(title = paste("Local Ancestry for", subj), 
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 20)))

ggsave(args$output, combined, width = 16, height = ceiling(length(plots)/4) * 4, dpi = 150)
message("Saved: ", args$output)
