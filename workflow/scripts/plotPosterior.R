library(argparse)
library(tidyverse) |> suppressPackageStartupMessages()

parser <- ArgumentParser(description = "Plot stacked posterior probabilities.")
parser$add_argument("--prob_file", type = "character", default = NULL,
    help = "Path to posterior_probabilities.tsv")
parser$add_argument("--out_dir", type = "character", default = NULL,
    help = "Output directory")
args <- parser$parse_args()

prob_df <- read_delim(args$prob_file, delim = "\t")

has_umap <- any(str_starts(colnames(prob_df), "umap_"))
has_vae <- any(str_starts(colnames(prob_df), "vae_"))
has_rfmix <- any(str_starts(colnames(prob_df), "rfmix_"))

available_models <- c("pca")
if (has_umap) available_models <- c(available_models, "umap")
if (has_vae) available_models <- c(available_models, "vae")
if (has_rfmix) available_models <- c(available_models, "rfmix")

for (model in available_models) {
    prob_cols <- prob_df |> select(IID, matches(paste0("^", model, "_"))) |> colnames()
    prob_cols <- setdiff(prob_cols, "IID")

    prob_model <- prob_df |>
        select(IID, all_of(prob_cols)) |>
        pivot_longer(cols = -IID, names_to = "ancestry", names_pattern = paste0(model, "_(.+)")) |>
        drop_na(value)

    avg_prob <- prob_model |>
        group_by(ancestry) |>
        summarise(avg_probability = mean(value, na.rm = TRUE), .groups = "drop") |>
        arrange(desc(avg_probability))

    sorted_ancestries <- avg_prob |> pull(ancestry)

    prob_model <- prob_model |>
        mutate(ancestry = factor(ancestry, levels = sorted_ancestries))

    subject_order <- prob_model |>
        select(IID, ancestry, value) |>
        pivot_wider(names_from = ancestry, values_from = value) |>
        arrange(!!!syms(rev(sorted_ancestries))) |>
        pull(IID)

    prob_model <- prob_model |>
        mutate(IID = factor(IID, levels = subject_order)) |>
        arrange(IID, ancestry)

    n_subjects <- length(unique(prob_model$IID))

    p <- prob_model |>
        ggplot(aes(x = IID, y = value, fill = ancestry)) +
        geom_col(position = "stack", width = 1) +
        theme_minimal() +
        theme(
            legend.position = "right",
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
        ) +
        xlab("Subjects") +
        ylab("Posterior Probability") +
        labs(fill = "Ancestry", title = paste0("Global Ancestry Proportions (", toupper(model), ")"))

    ggsave(file.path(args$out_dir, paste0("posterior_probability_stacked_", model, ".svg")),
        plot = p, width = 10, height = 7, units = "in")
}
