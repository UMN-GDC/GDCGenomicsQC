library(argparse)
library(tidyverse) |> suppressPackageStartupMessages()

parser <- ArgumentParser(description = "Plot ancestry classification space.")
parser$add_argument("--out_dir", type = "character", default = NULL,
    help = "Output directory containing classification files")
args <- parser$parse_args()

classification_df <- read_delim(file.path(args$out_dir, "ancestry_classifications.tsv"), delim = "\t")
sample_coords <- read_delim(file.path(args$out_dir, "sample_coords.tsv"), delim = "\t")
ref_data <- read_delim(file.path(args$out_dir, "ref_coords.tsv"), delim = "\t")

has_umap <- any(str_starts(colnames(sample_coords), "umap_"))
has_vae <- any(str_starts(colnames(sample_coords), "vae_"))
has_rfmix <- any(str_starts(colnames(classification_df), "rfmix_"))

sample_list <- list()
ref_list <- list()

sample_list$pca <- sample_coords |>
    select(IID, x = pc_1, y = pc_2) |>
    left_join(classification_df |> select(IID, pca_predicted, pca_confidence), by = "IID") |>
    rename(predicted = pca_predicted, confidence = pca_confidence) |>
    mutate(model = "pca")

ref_list$pca <- ref_data |>
    select(IID, x = pc_1, y = pc_2, POP) |>
    mutate(model = "pca")

if (has_umap && all(c("umap_1", "umap_2") %in% colnames(sample_coords))) {
    sample_list$umap <- sample_coords |>
        select(IID, x = umap_1, y = umap_2) |>
        left_join(classification_df |> select(IID, umap_predicted, umap_confidence), by = "IID") |>
        rename(predicted = umap_predicted, confidence = umap_confidence) |>
        mutate(model = "umap")

    ref_list$umap <- ref_data |>
        select(IID, x = umap_1, y = umap_2, POP) |>
        mutate(model = "umap")
} else {
    has_umap <- FALSE
}

if (has_vae && all(c("vae_mean1", "vae_mean2") %in% colnames(sample_coords))) {
    sample_list$vae <- sample_coords |>
        select(IID, x = vae_mean1, y = vae_mean2) |>
        left_join(classification_df |> select(IID, vae_predicted, vae_confidence), by = "IID") |>
        rename(predicted = vae_predicted, confidence = vae_confidence) |>
        mutate(model = "vae")

    ref_list$vae <- ref_data |>
        select(IID, x = vae_mean1, y = vae_mean2, POP) |>
        mutate(model = "vae")
} else {
    has_vae <- FALSE
}

available_models <- c("pca")
if (has_umap) available_models <- c(available_models, "umap")
if (has_vae) available_models <- c(available_models, "vae")
if (has_rfmix) available_models <- c(available_models, "rfmix")

rfmix_models <- available_models[available_models != "rfmix"]

sample_plot_df <- bind_rows(sample_list[rfmix_models])
ref_plot_df <- bind_rows(ref_list[rfmix_models])

if (nrow(sample_plot_df) > 0) {
    p <- ggplot() +
        stat_density_2d(data = ref_plot_df, aes(x = x, y = y, fill = POP),
            geom = "polygon", alpha = 0.25, contour = TRUE) +
        geom_point(data = sample_plot_df, aes(x = x, y = y, color = predicted),
            shape = 21, size = 2, stroke = 0.5) +
        facet_wrap(~model, ncol = 1, scales = "free") +
        theme_minimal() +
        theme(legend.position = "bottom")

    ggsave(file.path(args$out_dir, "ancestry_classification_space.svg"),
        plot = p, width = 1920, height = 800 * length(rfmix_models), units = "px")
}
