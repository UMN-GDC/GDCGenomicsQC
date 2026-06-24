#!/usr/bin/env Rscript

library(argparse)
suppressPackageStartupMessages(library(tidyverse))
library(mvtnorm)
library(jsonlite)

parser <- ArgumentParser(description = "Run N-ancestry phenotype simulation pipeline")
parser$add_argument("--anc-names", type = "character", required = TRUE,
    help = "Comma-separated ancestry names (e.g., AFR,EUR)")
parser$add_argument("--pgen-prefixes", type = "character", required = TRUE,
    help = "Comma-separated PGEN prefixes (no extension)")
parser$add_argument("--out-dirs", type = "character", required = TRUE,
    help = "Comma-separated per-ancestry output directories")
parser$add_argument("--corr-matrix", type = "character", required = TRUE,
    help = "NxN genetic correlation matrix as JSON, e.g. [[1.0,0.8],[0.8,1.0]]")
parser$add_argument("--n_sims", type = "integer", default = 10)
parser$add_argument("--seed", type = "integer", default = 42)
parser$add_argument("--heritability", type = "numeric", default = 0.4)
parser$add_argument("--maf", type = "numeric", default = 0.05)
parser$add_argument("--skip_thinning", type = "character", default = "true")
parser$add_argument("--thin_count_snps", type = "integer", default = 1000000)
parser$add_argument("--thin_count_inds", type = "integer", default = 10000)

args <- parser$parse_args()

anc_names <- strsplit(args$anc_names, ",")[[1]]
pgen_prefixes <- strsplit(args$pgen_prefixes, ",")[[1]]
out_dirs <- strsplit(args$out_dirs, ",")[[1]]
corr_matrix <- fromJSON(args$corr_matrix)

n_anc <- length(anc_names)
if (n_anc != nrow(corr_matrix) || n_anc != ncol(corr_matrix)) {
    stop(sprintf("corr_matrix dimension (%dx%d) does not match anc_names (n=%d)",
        nrow(corr_matrix), ncol(corr_matrix), n_anc))
}
if (n_anc != length(pgen_prefixes) || n_anc != length(out_dirs)) {
    stop("anc-names, pgen-prefixes, and out-dirs must have the same number of elements")
}

herit <- args$heritability
maf <- args$maf
n_sims <- args$n_sims
seed <- args$seed
skip_thinning <- tolower(args$skip_thinning) %in% c("true", "1")

cat("========================================\n")
cat("N-Ancestry Phenotype Simulation Pipeline\n")
cat("========================================\n")
cat(sprintf("Ancestries: %s\n", paste(anc_names, collapse = ", ")))
cat(sprintf("N = %d populations\n", n_anc))
cat(sprintf("Heritability: %.2f\n", herit))
cat(sprintf("MAF threshold: %.2f\n", maf))
cat(sprintf("Seed: %d\n", seed))
cat(sprintf("N sims: %d\n", n_sims))
cat(sprintf("Skip thinning: %s\n", skip_thinning))
cat("========================================\n\n")

for (d in out_dirs) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

cat("STEP 1: Thinning and preprocessing\n")
cat("-----------------------------------\n")

if (skip_thinning) {
    cat("Skipping thinning step\n")
    for (k in 1:n_anc) {
        wd <- out_dirs[k]
        system2("plink2", args = c("--pfile", pgen_prefixes[k],
            "--maf", as.character(maf),
            "--freq", "--out", file.path(wd, "sim_work")))
        if (k == 1) {
            system2("awk", args = c("{print $2}", paste0(pgen_prefixes[k], ".pvar")),
                stdout = file.path(wd, "snps.txt"))
        }
    }
} else {
    cat("Applying thinning\n")
    wd1 <- out_dirs[1]
    system2("plink2", args = c("--pfile", pgen_prefixes[1],
        "--thin-count", as.character(args$thin_count_snps),
        "--thin-indiv-count", as.character(args$thin_count_inds),
        "--maf", as.character(maf),
        "--freq", "--out", file.path(wd1, "sim_work")))
    system2("awk", args = c("{print $2}", paste0(pgen_prefixes[1], ".pvar")),
        stdout = file.path(wd1, "snps.txt"))

    for (k in 2:n_anc) {
        wd <- out_dirs[k]
        system2("plink2", args = c("--pfile", pgen_prefixes[k],
            "--extract", file.path(out_dirs[1], "snps.txt"),
            "--thin-indiv-count", as.character(args$thin_count_inds),
            "--maf", as.character(maf),
            "--freq", "--out", file.path(wd, "sim_work")))
    }
}

cat("\nSTEP 2: Creating base files\n")
cat("-----------------------------------\n")

for (k in 1:n_anc) {
    wd <- out_dirs[k]
    system2("plink2", args = c("--pfile", pgen_prefixes[k],
        "--keep-allele-order",
        "--make-bed", "--out", file.path(wd, "base")))

    fam_path <- file.path(wd, "base.fam")
    if (!file.exists(fam_path)) {
        stop(paste0("Base .fam not found: ", fam_path))
    }

    fam_data <- read_table(fam_path, col_names = FALSE)
    write_delim(fam_data[, 1:5], fam_path, col_names = FALSE, delim = " ")
    file.copy(fam_path, file.path(wd, "simulated.fam"), overwrite = TRUE)
    cat(sprintf("  %s: done\n", anc_names[k]))
}

cat("\nSTEP 3: Simulation\n")
cat("-----------------------------------\n")

simEffects <- function(out_dirs_list, anc_names_v, herit, corr_mat, maf) {
    n <- length(out_dirs_list)
    freq_list <- list()
    C_vals <- numeric(n)

    for (k in 1:n) {
        freq_path <- file.path(out_dirs_list[k], "sim_work.afreq")
        if (!file.exists(freq_path)) {
            freq_path <- file.path(out_dirs_list[k], "base.afreq")
        }
        df <- read_table(freq_path)
        C_vals[k] <- sum(df$ALT_FREQS > maf, na.rm = TRUE)
        colnames(df)[ncol(df)] <- paste0("ALT_FREQS_", k)
        freq_list[[k]] <- df
    }

    combined <- freq_list[[1]]
    if (n >= 2) {
        for (k in 2:n) {
            combined <- inner_join(combined, freq_list[[k]],
                by = c("#CHROM", "ID", "REF", "ALT"))
        }
    }

    for (k in 1:n) {
        fc <- paste0("ALT_FREQS_", k)
        combined[[paste0("sdEff_", k)]] <- case_when(
            is.na(combined[[fc]]) ~ 0,
            combined[[fc]] < maf ~ 0,
            .default = sqrt(herit / C_vals[k])
        )
    }

    cov_template <- matrix(0, n, n)
    for (k in 1:n) {
        for (l in 1:n) {
            if (k != l) {
                cov_template[k, l] <- corr_mat[k, l] * herit /
                    (sqrt(C_vals[k]) * sqrt(C_vals[l]))
            }
        }
    }

    betas <- matrix(NA, nrow = nrow(combined), ncol = n)
    for (i in 1:nrow(combined)) {
        cov_mat <- cov_template
        for (k in 1:n) {
            cov_mat[k, k] <- combined[[paste0("sdEff_", k)]][i]
        }
        if (any(is.na(cov_mat)) || any(is.nan(cov_mat))) {
            betas[i, ] <- rep(0, n)
        } else {
            betas[i, ] <- rmvnorm(n = 1, mean = rep(0, n), sigma = cov_mat)
        }
    }

    for (k in 1:n) {
        fc <- paste0("ALT_FREQS_", k)
        ec <- paste0("effs_", k)
        combined[[ec]] <- betas[, k] * sqrt(2 * combined[[fc]] * (1 - combined[[fc]]))
        esum <- sum(combined[[ec]]^2, na.rm = TRUE)
        if (esum > 0) {
            combined[[ec]] <- combined[[ec]] * sqrt(herit / esum)
        }
    }
    return(combined)
}

adjustPheno <- function(sscore_file, out_dir, herit) {
    pheno <- read_table(sscore_file)
    sd_score <- sd(pheno$SCORE1_AVG, na.rm = TRUE)
    effects <- read_delim(file.path(out_dir, "effects.csv"), delim = "\t")
    if (sd_score > 0) {
        effects[, 3] <- effects[, 3] / sd_score * sqrt(herit)
    }
    write_delim(effects, file.path(out_dir, "effects.csv"), delim = "\t")
}

finalizePheno <- function(sscore_file, out_dir, herit) {
    pheno <- read_table(sscore_file)
    pheno <- pheno %>%
        mutate(
            error = rnorm(n = nrow(pheno), 0, sqrt(1 - herit)),
            error = error / sd(error, na.rm = TRUE) * sqrt(1 - herit),
            pheno = SCORE1_AVG + error
        ) %>%
        select(-error)
    write_delim(pheno, file.path(out_dir, "simulated.pheno"), delim = "\t")
}

for (sim_i in 1:n_sims) {
    cat(sprintf("Running simulation %d of %d...\n", sim_i, n_sims))
    current_seed <- seed + sim_i - 1
    set.seed(current_seed)

    combined_effs <- simEffects(out_dirs, anc_names, herit, corr_matrix, maf)

    for (k in 1:n_anc) {
        wd <- out_dirs[k]
        eff_df <- combined_effs %>% select(ID, REF, !!sym(paste0("effs_", k)))
        colnames(eff_df)[3] <- "eff"
        write_delim(eff_df, file.path(wd, "effects.csv"), delim = "\t")

        system2("plink2", args = c("--pfile", pgen_prefixes[k],
            "--score", file.path(wd, "effects.csv"), "1", "2", "3", "header",
            "--out", file.path(wd, "sim_work")))

        adjustPheno(file.path(wd, "sim_work.sscore"), wd, herit)

        system2("plink2", args = c("--pfile", pgen_prefixes[k],
            "--score", file.path(wd, "effects.csv"), "1", "2", "3", "header",
            "--out", file.path(wd, "sim_work")))

        finalizePheno(file.path(wd, "sim_work.sscore"), wd, herit)

        pheno_data <- read_table(file.path(wd, "simulated.pheno"), col_names = TRUE)
        pheno_val <- pheno_data[[ncol(pheno_data)]]

        fam_path <- file.path(wd, "simulated.fam")
        fam_data <- read_table(fam_path, col_names = FALSE)
        fam_data <- fam_data %>% mutate(!!paste0("pheno_", sim_i) := pheno_val)
        write_delim(fam_data, fam_path, col_names = FALSE, delim = " ")
    }
    cat(sprintf("  Sim %d complete\n", sim_i))
}

cat("\nSTEP 4: Creating final BED files\n")
cat("-----------------------------------\n")

for (k in 1:n_anc) {
    wd <- out_dirs[k]
    sim_fam <- file.path(wd, "simulated.fam")
    tmp_fam <- paste0(sim_fam, ".tmp")
    file.copy(sim_fam, tmp_fam, overwrite = TRUE)

    system2("plink2", args = c("--pfile", pgen_prefixes[k],
        "--keep-allele-order",
        "--make-bed", "--out", file.path(wd, "simulated")))

    file.rename(tmp_fam, sim_fam)
    cat(sprintf("  %s: BED created\n", anc_names[k]))
}

cat("\n========================================\n")
cat("Simulation complete!\n")
for (k in 1:n_anc) {
    fam_final <- read_table(file.path(out_dirs[k], "simulated.fam"), col_names = FALSE)
    n_pheno <- ncol(fam_final) - 5
    cat(sprintf("  %s: %d samples, %d phenotype columns\n",
        anc_names[k], nrow(fam_final), n_pheno))
}
cat("========================================\n")