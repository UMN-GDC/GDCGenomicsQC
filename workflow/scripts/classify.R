library(argparse)
library(tidyverse) |> suppressPackageStartupMessages()

parser <- ArgumentParser(description = "Classify ancestry based on posterior probabilities.")
parser$add_argument("--out", type = "character", default = NULL,
    help = "Output directory")
parser$add_argument("--threshold", type = "double", default = 0.8,
    help = "Threshold for ancestry classification (default: 0.8)")
parser$add_argument("--model", type = "character", default = "pca",
    help = "Model to use for classification (default: pca)")
args <- parser$parse_args()

prob_df <- read_delim(file.path(args$out, "posterior_probabilities.tsv"), delim = "\t")

has_umap <- any(colnames(prob_df) %>% str_starts("umap_"))
has_vae <- any(colnames(prob_df) %>% str_starts("vae_"))
has_rfmix <- any(colnames(prob_df) %>% str_starts("rfmix_"))

result <- prob_df |> select(IID)

for (model in c("pca", "umap", "vae", "rfmix")) {
    prob_cols <- prob_df |> select(starts_with(paste0(model, "_"))) |> colnames()
    if (length(prob_cols) == 0) next

    prob_matrix <- prob_df |> select(all_of(prob_cols)) |> as.matrix()
    max_idx <- max.col(prob_matrix, ties.method = "first")
    ancestry_names <- gsub(paste0(model, "_"), "", prob_cols)

    result[[paste0(model, "_predicted")]] <- ancestry_names[max_idx]
    result[[paste0(model, "_confidence")]] <- apply(prob_matrix, 1, max)
    result[[paste0(model, "_predicted")]] <- ifelse(
        result[[paste0(model, "_confidence")]] >= args$threshold,
        result[[paste0(model, "_predicted")]],
        "uncertain"
    )
}

write_delim(result, file.path(args$out, "ancestry_classifications.tsv"), delim = "\t")

model <- args$model
predicted_col <- paste0(model, "_predicted")
confidence_col <- paste0(model, "_confidence")

unique_ancestries <- result |>
    filter(.data[[confidence_col]] >= args$threshold) |>
    pull(predicted_col) |>
    unique()

for (anc in unique_ancestries) {
    result |>
        filter(.data[[predicted_col]] == anc, .data[[confidence_col]] >= args$threshold) |>
        select(IID) |>
        write_delim(file.path(args$out, paste0("keep_", anc, ".txt")), delim = "\t", col_names = FALSE)
}

result |>
    filter(.data[[confidence_col]] < args$threshold | .data[[predicted_col]] == "uncertain" | is.na(.data[[confidence_col]])) |>
    select(IID) |>
    write_delim(file.path(args$out, "keep_Other.txt"), delim = "\t", col_names = FALSE)

classification_df <- result
sample_coords <- read_delim(file.path(args$out, "sample_coords.tsv"), delim = "\t")
ref_data <- read_delim(file.path(args$out, "ref_coords.tsv"), delim = "\t")

sample_long <- list()
ref_long <- list()

sample_long$pca <- sample_coords |>
    select(IID, x = pc_1, y = pc_2) |>
    left_join(classification_df |> select(IID, pca_predicted, pca_confidence),
        by = "IID") |>
    rename(predicted = pca_predicted, confidence = pca_confidence) |>
    mutate(model = "pca")

ref_long$pca <- ref_data |>
    select(IID, x = pc_1, y = pc_2, POP) |>
    mutate(model = "pca")

if (has_umap & all(c("umap_1", "umap_2") %in% colnames(sample_coords))) {
    sample_long$umap <- sample_coords |>
        select(IID, x = umap_1, y = umap_2) |>
        left_join(classification_df |> select(IID, umap_predicted, umap_confidence),
            by = "IID") |>
        rename(predicted = umap_predicted, confidence = umap_confidence) |>
        mutate(model = "umap")

    ref_long$umap <- ref_data |>
        select(IID, x = umap_1, y = umap_2, POP) |>
        mutate(model = "umap")
} else {
    has_umap <- FALSE
}

if (has_vae & all(c("vae_mean1", "vae_mean2") %in% colnames(sample_coords))) {
    sample_long$vae <- sample_coords |>
        select(IID, x = vae_mean1, y = vae_mean2) |>
        left_join(classification_df |> select(IID, vae_predicted, vae_confidence),
            by = "IID") |>
        rename(predicted = vae_predicted, confidence = vae_confidence) |>
        mutate(model = "vae")

    ref_long$vae <- ref_data |>
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

sample_plot_df <- bind_rows(sample_long[rfmix_models])
ref_plot_df <- bind_rows(ref_long[rfmix_models])

if (nrow(sample_plot_df) > 0) {
    p <- ggplot() +
        stat_density_2d(data = ref_plot_df, aes(x = x, y = y, fill = POP),
            geom = "polygon", alpha = 0.25, contour = TRUE) +
        geom_point(data = sample_plot_df, aes(x = x, y = y, color = predicted),
            shape = 21, size = 2, stroke = 0.5) +
        facet_wrap(~model, ncol = 1, scales = "free") +
        theme_minimal() +
        theme(legend.position = "bottom")

    ggsave(file.path(args$out, "ancestry_classification_space.svg"),
        plot = p, dpi = 300, width = 9, height = 5 * length(rfmix_models))
}
