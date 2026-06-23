library(argparse)
library(tidyverse) |> suppressPackageStartupMessages()
library(randomForest)

parser <- ArgumentParser(description = "Train ancestry RF models and predict probabilities.")
parser$add_argument("--eigen_ref", type = "character", default = NULL,
    help = "Filepath to the eigenvector file (plink style)")
parser$add_argument("--eigen_sample", type = "character", default = NULL,
    help = "Filepath to the eigenvector file (plink style)")
parser$add_argument("--labels", type = "character", default = NULL,
    help = "Path to labels file. Assumes file columns, FID, IID, Population")
parser$add_argument("--vae", type = "character", default = NULL,
    help = "Filepath to the vae file (popvae style)")
parser$add_argument("--umap_ref", type = "character", default = NULL,
    help = "Filepath to the umap file")
parser$add_argument("--umap_sample", type = "character", default = NULL,
    help = "Filepath to the umap file")
parser$add_argument("--rfmix_global", type = "character", default = NULL,
    help = "Filepath to global RFMix ancestry output (ancestry_full.txt)")
parser$add_argument("--out", type = "character", default = NULL,
    help = "Output directory")
parser$add_argument("--rseed",
    type = "integer", default = as.integer(Sys.time()),
    help = "Specify the desired seed. Default system time")
args <- parser$parse_args()

set.seed(args$rseed)

read_vae_coords <- function(path) {
    vae_df <- read_table(path, col_names = TRUE, show_col_types = FALSE)
    names(vae_df) <- make.names(names(vae_df), unique = TRUE)

    iid_candidates <- c("IID", "SampleID", "sample_id", "sample", "id", "ID")
    iid_col <- iid_candidates[iid_candidates %in% names(vae_df)][1]
    if (is.na(iid_col)) {
        iid_col <- names(vae_df)[ncol(vae_df)]
    }

    names(vae_df)[names(vae_df) == iid_col] <- "IID"
    names(vae_df) <- ifelse(names(vae_df) == "IID", "IID", paste0("vae_", names(vae_df)))

    if (!("vae_mean1" %in% names(vae_df)) || !("vae_mean2" %in% names(vae_df))) {
        numeric_cols <- names(vae_df)[vapply(vae_df, is.numeric, logical(1))]
        numeric_cols <- setdiff(numeric_cols, "IID")
        if (length(numeric_cols) < 2) {
            stop("VAE coordinate file must contain IID plus at least two numeric latent-coordinate columns.")
        }
        vae_df$vae_mean1 <- vae_df[[numeric_cols[1]]]
        vae_df$vae_mean2 <- vae_df[[numeric_cols[2]]]
    }

    vae_df |>
        select(IID, starts_with("vae_")) |>
        distinct(IID, .keep_all = TRUE)
}

fit_and_predict_ancestry_models <- function(
    ref_labels,
    eigen_ref,
    eigen_sample,
    umap_ref = NULL,
    umap_sample = NULL,
    vae_ref = NULL,
    rfmix_global = NULL,
    out_dir
) {
    ref <- read_table(ref_labels) |>
        select(FID = FamilyID, IID = SampleID, POP = Superpopulation)

    ancestries <- unique(ref$POP)

    PCs <- read_table(eigen_ref, col_names = TRUE, show_col_types = FALSE)

    for (col in c("ALLELE_CT", "NAMED_ALLELE_DOSAGE_SUM", "FID", "#FID")) {
        if (col %in% colnames(PCs)) {
            PCs <- PCs |> select(-all_of(col))
        }
    }

    iid_col <- if ("#IID" %in% colnames(PCs)) "#IID" else "IID"
    names(PCs)[names(PCs) == iid_col] <- "IID"
    colnames(PCs)[-1] <- paste0("pc_", 1:(ncol(PCs) - 1))
    ref <- full_join(ref, PCs, by = c("IID")) |> drop_na(pc_1)

    pcMod <- randomForest::randomForest(
        formula = factor(POP) ~ pc_1 + pc_2 + pc_3 + pc_4 + pc_5 + pc_6 + pc_7 + pc_8 + pc_9 + pc_10,
        data = ref
    )
    saveRDS(pcMod, file.path(out_dir, "RFpc.Rds"))

    sampleDF <- read_table(eigen_sample, col_names = TRUE, show_col_types = FALSE)

    for (col in c("ALLELE_CT", "NAMED_ALLELE_DOSAGE_SUM", "FID", "#FID")) {
        if (col %in% colnames(sampleDF)) {
            sampleDF <- sampleDF |> select(-all_of(col))
        }
    }

    iid_col <- if ("#IID" %in% colnames(sampleDF)) "#IID" else "IID"
    names(sampleDF)[names(sampleDF) == iid_col] <- "IID"
    colnames(sampleDF)[-1] <- paste0("pc_", 1:(ncol(sampleDF) - 1))

    pc_probs <- randomForest:::predict.randomForest(pcMod, sampleDF, type = "prob")
    result_df <- sampleDF |> select(IID) |> as_tibble()
    sample_coords_df <- sampleDF

    for (anc in ancestries) {
        result_df[[paste0("pca_", anc)]] <- pc_probs[, anc]
    }

    has_umap <- FALSE
    if (!is.null(umap_ref)) {
        has_umap <- TRUE
        umap_ref_df <- read_csv(umap_ref)
        colnames(umap_ref_df) <- c("IID", str_replace(colnames(umap_ref_df)[-c(1)], "UMAP", "umap_"))
        ref <- full_join(ref, umap_ref_df, by = c("IID")) |> drop_na(umap_1)

        umapMod <- randomForest::randomForest(
            formula = factor(POP) ~ umap_1 + umap_2,
            data = ref
        )
        saveRDS(umapMod, file.path(out_dir, "RFumap.Rds"))

        umap_sample_df <- read_csv(umap_sample)
        colnames(umap_sample_df) <- c("IID", str_replace(colnames(umap_sample_df)[-c(1)], "UMAP", "umap_"))

        sampleDF_umap <- sampleDF |>
            inner_join(umap_sample_df, by = "IID")
        sample_coords_df <- sample_coords_df |>
            left_join(umap_sample_df, by = "IID")

        umap_probs <- randomForest:::predict.randomForest(umapMod, sampleDF_umap, type = "prob")

        umap_result <- sampleDF_umap |>
            select(IID) |>
            as_tibble()
        for (anc in ancestries) {
            umap_result[[paste0("umap_", anc)]] <- umap_probs[, anc]
        }

        result_df <- result_df |>
            left_join(umap_result, by = "IID")
    }

    has_vae <- FALSE
    if (!is.null(vae_ref)) {
        has_vae <- TRUE
        vae_ref_df <- read_vae_coords(vae_ref)
        ref <- full_join(ref, vae_ref_df, by = c("IID")) |> drop_na(vae_mean1, vae_mean2)

        vaeMod <- randomForest::randomForest(
            formula = factor(POP) ~ vae_mean1 + vae_mean2,
            data = ref
        )
        saveRDS(vaeMod, file.path(out_dir, "RFvae.Rds"))

        vae_sample_df <- vae_ref_df

        sampleDF_vae <- sampleDF |>
            inner_join(vae_sample_df, by = "IID")
        sample_coords_df <- sample_coords_df |>
            left_join(vae_sample_df, by = "IID")

        if (nrow(sampleDF_vae) == 0) {
            warning("VAE coordinates not found in sample data. Skipping VAE prediction.")
            has_vae <- FALSE
        } else {
            vae_probs <- randomForest:::predict.randomForest(vaeMod, sampleDF_vae, type = "prob")

            vae_result <- sampleDF_vae |>
                select(IID) |>
                as_tibble()
            for (anc in ancestries) {
                vae_result[[paste0("vae_", anc)]] <- vae_probs[, anc]
            }

            result_df <- result_df |>
                left_join(vae_result, by = "IID")
        }
    }

    has_rfmix <- FALSE
    if (!is.null(rfmix_global)) {
        has_rfmix <- TRUE
        rfmix_df <- read.table(rfmix_global, header = TRUE, check.names = FALSE) |>
            as_tibble()

        rfmix_ancestries <- colnames(rfmix_df) |> setdiff("IID")

        result_df <- result_df |>
            left_join(rfmix_df, by = "IID")

        for (anc in rfmix_ancestries) {
            if (!(paste0("rfmix_", anc) %in% colnames(result_df))) {
                result_df[[paste0("rfmix_", anc)]] <- result_df[[anc]]
                result_df[[anc]] <- NULL
            }
        }
    }

    list(
        probabilities = result_df,
        sample_coords = sample_coords_df,
        ref_data = ref,
        ancestries = ancestries,
        has_umap = has_umap,
        has_vae = has_vae,
        has_rfmix = has_rfmix
    )
}

prob_results <- fit_and_predict_ancestry_models(
    ref_labels = args$labels,
    eigen_ref = args$eigen_ref,
    eigen_sample = args$eigen_sample,
    umap_ref = args$umap_ref,
    umap_sample = args$umap_sample,
    vae_ref = args$vae,
    rfmix_global = args$rfmix_global,
    out_dir = args$out
)

prob_results$probabilities |>
    relocate(IID) |>
    write_delim(file.path(args$out, "classificationProbabilities.tsv"), delim = "\t")

prob_results$sample_coords |>
    relocate(IID) |>
    write_delim(file.path(args$out, "sample_coords.tsv"), delim = "\t")

ref_for_plot <- prob_results$ref_data |>
    select(IID, POP, starts_with("pc_"), starts_with("umap_"), starts_with("vae_")) |>
    relocate(IID)
ref_for_plot |>
    write_delim(file.path(args$out, "ref_coords.tsv"), delim = "\t")
