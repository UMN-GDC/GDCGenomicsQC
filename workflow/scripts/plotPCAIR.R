library(argparse)
library(tidyverse) |> suppressPackageStartupMessages()

parser <- ArgumentParser(description = "Plot PC-AiR internal PCA results.")
parser$add_argument("--coords", type = "character", default = NULL,
    help = "Path to pcair_coordinates.tsv")
parser$add_argument("--out", type = "character", default = NULL,
    help = "Output SVG path")
parser$add_argument("--color-col", type = "character", default = "None",
    help = "Column name in phenotype file to color by")
parser$add_argument("--pheno-file", type = "character", default = "None",
    help = "Path to phenotype file for coloring")
args <- parser$parse_args()

coords <- read_delim(args$coords, delim = "\t")

if (args$color_col != "None" && args$pheno_file != "None" && file.exists(args$pheno_file)) {
    pheno <- read_tsv(args$pheno_file, show_col_types = FALSE)
    if (args$color_col %in% colnames(pheno)) {
        coords <- coords |> left_join(pheno |> select(IID, all_of(args$color_col)), by = "IID")
        color_var <- args$color_col
    } else {
        coords$color_group <- "all"
        color_var <- "color_group"
    }
} else {
    coords$color_group <- "all"
    color_var <- "color_group"
}

p <- ggplot(coords, aes(x = PC1, y = PC2, color = .data[[color_var]])) +
    geom_point(alpha = 0.6, size = 2) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(title = "PC-AiR: Internal Sample PCs")

ggsave(args$out, plot = p, dpi = 300, width = 8, height = 6)
