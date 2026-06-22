library(argparse)
library(tidyverse) |> suppressPackageStartupMessages()

parser <- ArgumentParser(description = "Plot classification probability density by ancestry.")
parser$add_argument("--prob_file", type = "character", default = NULL,
    help = "Path to classificationProbabilities.tsv")
parser$add_argument("--out_dir", type = "character", default = NULL,
    help = "Output directory")
parser$add_argument("--model", type = "character", default = "pca",
    help = "Model prefix to plot (default: pca)")
args <- parser$parse_args()

prob_df <- read_delim(args$prob_file, delim = "\t")

model <- args$model
prob_cols <- prob_df |> select(IID, starts_with(paste0(model, "_"))) |> colnames()
prob_cols <- setdiff(prob_cols, "IID")

if (length(prob_cols) == 0) {
    stop("No columns found for model prefix: ", model)
}

df <- prob_df |>
    select(IID, all_of(prob_cols)) |>
    pivot_longer(cols = -IID, names_to = "ANC", values_to = "proportion") |>
    mutate(ANC = str_remove(ANC, paste0("^", model, "_")))

p <- df |>
    ggplot(aes(x = proportion, col = ANC)) +
    geom_density(alpha = 0.3) +
    scale_x_continuous(limits = c(0.5, 1)) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(
        x = "Classification Probability",
        y = "Density",
        col = "Ancestry",
        title = paste0("Classification Probability Density (", toupper(model), ")")
    )

ggsave(file.path(args$out_dir, paste0("classificationProbability_density_", model, ".svg")),
    plot = p, width = 10, height = 6, units = "in")
