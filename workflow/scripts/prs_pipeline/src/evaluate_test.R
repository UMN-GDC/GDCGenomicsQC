#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop("Usage: evaluate_test.R <profile> <pheno> <pca> <results_out> <scores_out>")
}

profile_file <- args[1]
pheno_file <- args[2]
pca_file <- args[3]
results_out <- args[4]
scores_out <- args[5]

# Helper: detect header (col 2 of first row has letters)
detect_header <- function(file) {
  first <- read.table(file, nrows = 1, stringsAsFactors = FALSE)
  grepl("[A-Za-z]", first[1, 2])
}

# Read profile (plink --score output with header: FID IID PHENO CNT CNT2 SCORE)
profile <- read.table(profile_file, header = TRUE)

# Read phenotype
pheno_has_header <- detect_header(pheno_file)
phenotype <- read.table(pheno_file, header = pheno_has_header)
colnames(phenotype)[1:3] <- c("FID", "IID", "phenotype")

# Read PCs
pca_has_header <- detect_header(pca_file)
pcs <- read.table(pca_file, header = pca_has_header)
colnames(pcs)[1:2] <- c("FID", "IID")
if (ncol(pcs) > 2) {
  colnames(pcs)[3:ncol(pcs)] <- paste0("PC", 1:(ncol(pcs) - 2))
}

# Merge
data <- merge(phenotype, pcs, by = c("FID", "IID"))
data <- merge(data, profile[, c("FID", "IID", "SCORE")], by = c("FID", "IID"))

if (nrow(data) == 0) {
  stop("Merge produced 0 rows — FID/IID mismatch between profile, phenotype, and PCA files")
}

# Null model (PCs only)
null_vars <- grep("^PC", colnames(data), value = TRUE)
null_formula <- as.formula(paste("phenotype ~", paste(null_vars, collapse = " + ")))
null_model <- lm(null_formula, data = data)
null_r2 <- summary(null_model)$r.squared

# Full model (PCs + SCORE)
full_formula <- as.formula(paste("phenotype ~ SCORE +", paste(null_vars, collapse = " + ")))
full_model <- lm(full_formula, data = data)
full_r2 <- summary(full_model)$r.squared
prs_r2 <- full_r2 - null_r2

prs_coef <- summary(full_model)$coeff["SCORE", , drop = FALSE]
prs_beta <- prs_coef[1, "Estimate"]
prs_se   <- prs_coef[1, "Std. Error"]
prs_p    <- prs_coef[1, "Pr(>|t|)"]

# Write results
results <- data.frame(R2 = prs_r2, P = prs_p, BETA = prs_beta, SE = prs_se)
write.table(results, file = results_out, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Write individual scores
scores <- data.frame(FID = data$FID, IID = data$IID, SCORE = data$SCORE)
write.table(scores, file = scores_out, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

cat("Test evaluation complete: R2 =", round(prs_r2, 5), "N =", nrow(data), "\n")
