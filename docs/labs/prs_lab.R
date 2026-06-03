# Lab: Polygenic Risk Score Analysis
#
# Compare PRS methods (C+T, LDpred2, lassosum2, PRSice-2) and
# calculate incremental R2 using a null + full model approach.
#
# Prerequisites:
# - Completed PRS pipeline run with output in output_path
# - R (>= 4.0) with tidyverse and corrplot installed

conda_env <- '/scratch.global/R25_files/r25_r_env/lib/R/library'
.libPaths(unique(c(.Library.site, .Library)))
.libPaths(c(conda_env, .libPaths()))

library(tidyverse)
library(data.table)

output_path <- "/path/to/your/results"

# 1. Load C+T scores
ct_path <- file.path(output_path, "prs_pipeline/CT/temp/temp.0.X.profile")
ct_data <- fread(ct_path)

# 2. Load LDpred2 scores
ldp_path <- file.path(output_path, "prs_pipeline/LDpred2/prs_method_individual_scores.txt")
ldp_data <- fread(ldp_path)

# 3. Load lassosum2 scores
lasso_path <- file.path(output_path, "prs_pipeline/lassosum2/prs_method_grid_params.csv")
lasso_data <- fread(lasso_path)

# Combine scores
combined_scores <- data.frame(
  IID = ct_data$IID,
  CT = ct_data$SCORE,
  LDpred2 = ldp_data$PRS_grid,
  lassosum2 = lasso_data$score
)

# Correlation matrix
cor_matrix <- cor(combined_scores[, -1], use = "complete.obs")
print(cor_matrix)

library(corrplot)
corrplot(cor_matrix, method = "color", addCoef.col = "black")

# Load phenotype data
pheno_path <- "/path/to/your/phenotypes.txt"
pheno_data <- fread(pheno_path)

final_data <- merge(combined_scores, pheno_data, by = "IID")
final_data$Pheno <- as.numeric(final_data$Pheno)

# Calculate incremental R2
covariates <- "Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5"

calculate_incremental_r2 <- function(prs_col, data) {
  null_formula <- as.formula(paste("Pheno ~", covariates))
  null_model <- lm(null_formula, data = data)
  r2_null <- summary(null_model)$r.squared

  full_formula <- as.formula(paste("Pheno ~", covariates, "+", prs_col))
  full_model <- lm(full_formula, data = data)
  r2_full <- summary(full_model)$r.squared

  r2_full - r2_null
}

r2_results <- data.frame(
  Method = c("CT", "LDpred2", "lassosum2"),
  Incremental_R2 = c(
    calculate_incremental_r2("CT", final_data),
    calculate_incremental_r2("LDpred2", final_data),
    calculate_incremental_r2("lassosum2", final_data)
  )
)
print(r2_results)

# Visualize performance
ggplot(r2_results, aes(x = reorder(Method, -Incremental_R2), y = Incremental_R2, fill = Method)) +
  geom_bar(stat = "identity", color = "black", show.legend = FALSE) +
  theme_minimal() +
  labs(
    title = "PRS Performance Comparison",
    subtitle = "Measured by Incremental R-squared",
    x = "Statistical Method",
    y = expression("Incremental R"^2)
  ) +
  scale_fill_brewer(palette = "Set2")

# Scatter plot for top method (LDpred2)
ggplot(final_data, aes(x = LDpred2, y = Pheno)) +
  geom_point(alpha = 0.3, color = "steelblue") +
  geom_smooth(method = "lm", color = "darkred") +
  theme_light() +
  labs(
    title = "Observed Phenotype vs. LDpred2 Score",
    x = "Polygenic Risk Score (LDpred2)",
    y = "Observed Phenotype Value"
  )
