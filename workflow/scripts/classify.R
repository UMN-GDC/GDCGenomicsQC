library(argparse)
library(tidyverse) |> suppressPackageStartupMessages()

parser <- ArgumentParser(description = "Classify ancestry based on classification probabilities.")
parser$add_argument("--out", type = "character", default = NULL,
    help = "Output directory")
parser$add_argument("--threshold", type = "double", default = 0.8,
    help = "Threshold for ancestry classification (default: 0.8)")
parser$add_argument("--model", type = "character", default = "pca",
    help = "Model to use for classification (default: pca)")
args <- parser$parse_args()

prob_df <- read_delim(file.path(args$out, "classificationProbabilities.tsv"), delim = "\t")

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

prob_cols <- prob_df |> select(starts_with(paste0(model, "_"))) |> colnames()
all_ancestries <- gsub(paste0(model, "_"), "", prob_cols)

for (anc in all_ancestries) {
    keep_df <- result |>
        filter(.data[[predicted_col]] == anc, .data[[confidence_col]] >= args$threshold) |>
        select(IID, any_of("FID"))
    if (!"FID" %in% colnames(keep_df)) {
        keep_df <- keep_df |> mutate(FID = IID)
    }
    keep_df |> relocate(FID) |> write_delim(
        file.path(args$out, paste0("keep_", anc, ".txt")), delim = "\t", col_names = FALSE
    )
}

other_df <- result |>
    filter(.data[[confidence_col]] < args$threshold | .data[[predicted_col]] == "uncertain" | is.na(.data[[confidence_col]])) |>
    select(IID, any_of("FID"))
if (!"FID" %in% colnames(other_df)) {
    other_df <- other_df |> mutate(FID = IID)
}
other_df |> relocate(FID) |>
    write_delim(file.path(args$out, "keep_Other.txt"), delim = "\t", col_names = FALSE)


