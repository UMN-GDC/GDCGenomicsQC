library(argparse)
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(magrittr)
library(ranger)

parser <- ArgumentParser(description = "Plot ancestry classification space.")
parser$add_argument("--out_dir", type = "character", default = NULL,
    help = "Output directory containing classification files")
parser$add_argument("--threshold", type = "double", default = 0.8,
    help = "Probability threshold for classification (default: 0.8)")
parser$add_argument("--model", type = "character", default = "pca",
    help = "Model prefix (default: pca)")
parser$add_argument("--rf_model", type = "character", default = NULL,
    help = "Path to RF model .Rds file")
args <- parser$parse_args()

classification_df <- read_delim(file.path(args$out_dir, "ancestry_classifications.tsv"), delim = "\t")
sample_coords <- read_delim(file.path(args$out_dir, "sample_coords.tsv"), delim = "\t") |>
    mutate(across(any_of(c("pc_1", "pc_2", "umap_1", "umap_2", "vae_mean1", "vae_mean2")), as.numeric))
ref_data <- read_delim(file.path(args$out_dir, "ref_coords.tsv"), delim = "\t") |>
    mutate(across(any_of(c("pc_1", "pc_2", "umap_1", "umap_2", "vae_mean1", "vae_mean2")), as.numeric))

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

contour_list <- list()

if (!is.null(args$rf_model) && file.exists(args$rf_model) && "pca" %in% rfmix_models) {
    rf_model <- readRDS(args$rf_model)
    rf_vars <- all.vars(formula(rf_model))

    pca_ref <- ref_plot_df |> filter(model == "pca")

    x_range <- range(pca_ref$x, na.rm = TRUE)
    y_range <- range(pca_ref$y, na.rm = TRUE)

    grid <- expand.grid(
        pc_1 = seq(x_range[1], x_range[2], length.out = 150),
        pc_2 = seq(y_range[1], y_range[2], length.out = 150)
    )

    for (var in setdiff(rf_vars, c("pc_1", "pc_2"))) {
        if (var %in% colnames(ref_data)) {
            grid[[var]] <- median(ref_data[[var]], na.rm = TRUE)
        }
    }

    probs <- predict(rf_model, grid)$predictions
    grid$max_prob <- apply(probs, 1, max)

    contour_list$pca <- grid |>
        rename(x = pc_1, y = pc_2) |>
        mutate(model = "pca")
}

if (has_umap && "umap" %in% rfmix_models) {
    umap_rf_path <- file.path(args$out_dir, "RFumap.Rds")
    if (file.exists(umap_rf_path)) {
        umap_model <- readRDS(umap_rf_path)
        umap_ref <- ref_plot_df |> filter(model == "umap")

        x_range <- range(umap_ref$x, na.rm = TRUE)
        y_range <- range(umap_ref$y, na.rm = TRUE)

        grid <- expand.grid(
            umap_1 = seq(x_range[1], x_range[2], length.out = 150),
            umap_2 = seq(y_range[1], y_range[2], length.out = 150)
        )

        probs <- predict(umap_model, grid)$predictions
        grid$max_prob <- apply(probs, 1, max)

        contour_list$umap <- grid |>
            rename(x = umap_1, y = umap_2) |>
            mutate(model = "umap")
    }
}

contour_df <- bind_rows(contour_list)

if (nrow(sample_plot_df) > 0) {
    p <- ggplot() +
        geom_density_2d(data = ref_plot_df, aes(x = x, y = y, color = POP),
            breaks = c(0.7)) +
        geom_point(data = sample_plot_df, aes(x = x, y = y, color = predicted),
            alpha = 0.3, size = 0.5)

    if (nrow(contour_df) > 0) {
        p <- p + geom_contour(data = contour_df, aes(x = x, y = y, z = max_prob),
            breaks = args$threshold, color = "red", linewidth = 0.6)
    }

    p <- p +
        facet_wrap(~model, ncol = 1, scales = "free") +
        theme_minimal() +
        theme(legend.position = "bottom")

    ggsave(file.path(args$out_dir, "ancestry_classification_space.svg"),
        plot = p, width = 1920, height = 800 * length(rfmix_models), units = "px")
}
