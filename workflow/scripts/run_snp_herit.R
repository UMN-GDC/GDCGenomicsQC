#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 9) {
    stop("Usage: Rscript run_snp_herit.R <out_dir> <pcrelate_file> <pcaobj_file> <unrelated_ids_file> <pheno_file> <covar_file> <out_file> <npc> <mpheno> <method>\n")
}

out_dir <- args[1]
pcrelate_file <- args[2]
pcaobj_file <- args[3]
unrels_file <- args[4]
pheno_file <- args[5]
covar_file <- args[6]
out_file <- args[7]
npc <- as.integer(args[8])
mpheno <- as.integer(args[9])
method <- args[10]

suppressPackageStartupMessages({
    library(GENESIS)
    library(SeqArray)
    library(SeqVarTools)
    library(dplyr)
    library(purrr)
})

cat("Loading PC-Relate kinship matrix...\n")
pcrelate <- readRDS(pcrelate_file)

cat("Loading PC-AiR PCA object...\n")
pcaobj <- readRDS(pcaobj_file)

cat("Reading unrelated sample IDs...\n")
unrels <- read.table(unrels_file, stringsAsFactors = FALSE)$V1
cat("Found", length(unrels), "unrelated samples\n")

pheno_ext <- tools::file_ext(pheno_file)
if (pheno_ext == "fam") {
    cat("Reading phenotype file as FAM format...\n")
    pheno <- read.table(pheno_file, header = FALSE, stringsAsFactors = FALSE)
    colnames(pheno) <- c("FID", "IID", "PID", "MID", "Sex", "Pheno")
    pheno <- pheno[, c("IID", "Pheno")]
    colnames(pheno) <- c("IID", "phenotype")
} else {
    cat("Reading phenotype file as CSV...\n")
    pheno <- read.csv(pheno_file, stringsAsFactors = FALSE)
}

cat("Phenotype columns:", colnames(pheno), "\n")

cat("Reading covariate file(s)...\n")
covar_files <- unlist(strsplit(covar_file, ","))
cat("Covariate file(s):", covar_files, "\n")

covar_list <- lapply(covar_files, function(f) {
    cat("Reading:", f, "\n")
    read.csv(f, stringsAsFactors = FALSE)
})

if (length(covar_list) > 1) {
    cat("Merging", length(covar_list), "covariate files on IID...\n")
    covar <- Reduce(function(x, y) merge(x, y, by = "IID", all = FALSE, suffixes = c("", ".y")), covar_list)
    covar <- covar[, !grepl("\\.y$", colnames(covar))]
} else {
    covar <- covar_list[[1]]
}

cat("Covariate columns:", colnames(covar), "\n")

cat("Matching samples...\n")
pheno <- pheno %>% filter(IID %in% unrels)
covar <- covar %>% filter(IID %in% unrels)

common_ids <- intersect(pheno$IID, covar$IID)
pheno <- pheno %>% filter(IID %in% common_ids)
covar <- covar %>% filter(IID %in% common_ids)

pheno <- pheno %>% arrange(match(IID, common_ids))
covar <- covar %>% arrange(match(IID, common_ids))

cat("Final sample count:", length(common_ids), "\n")

cat("Computing GRM from PC-Relate...\n")
grm <- pcrelateToMatrix(pcrelate, pcaobj = pcaobj, shrink = TRUE)
dimnames(grm) <- list(rownames(grm), colnames(grm))

pc_matrix <- pcaobj$vectors[, 1:npc, drop = FALSE]
rownames(pc_matrix) <- pcaobj$sample.id
pc_matrix <- pc_matrix[common_ids, , drop = FALSE]
grm <- grm[common_ids, common_ids]
grm <- grm[match(common_ids, rownames(grm)), match(common_ids, colnames(grm))]

cat("Running heritability estimation with method:", method, "\n")

if (method == "AdjHE") {
    cat("Using Haseman-Elston regression...\n")
    
    y <- pheno$phenotype
    names(y) <- pheno$IID
    
    X <- cbind(rep(1, length(common_ids)), pc_matrix)
    
    n <- length(y)
    p <- ncol(X)
    
    y_centered <- y - mean(y)
    X_centered <- scale(X, center = TRUE, scale = FALSE)
    
    resid <- lm(y ~ X_centered - 1)$residuals
    
    cat("Computing covariance matrix...\n")
    Sigma <- grm
    diag(Sigma) <- diag(Sigma) + 0.01
    
    cat("Fitting mixed model with EM...\n")
    sigma_g <- var(resid) * 0.5
    sigma_e <- var(resid) * 0.5
    
    for (iter in 1:50) {
        V <- sigma_g * Sigma + sigma_e * diag(n)
        V_inv <- solve(V)
        
        XtVinvX <- t(X_centered) %*% V_inv %*% X_centered
        XtVinvy <- t(X_centered) %*% V_inv %*% resid
        
        beta <- solve(XtVinvX, XtVinvy)
        
        resid <- resid - X_centered %*% beta
        
        tr_sigma <- sum(diag(Sigma %*% V_inv))
        tr_sigma2 <- sum(diag((Sigma %*% V_inv) %*% (Sigma %*% V_inv)))
        
        sigma_g_new <- (t(resid) %*% V_inv %*% Sigma %*% V_inv %*% resid / n) * 
                       (tr_sigma / tr_sigma2)
        sigma_e_new <- (t(resid) %*% V_inv %*% resid / n) * 
                       ((n - p) / tr_sigma2)
        
        sigma_g_new <- max(sigma_g_new, 0.001)
        sigma_e_new <- max(sigma_e_new, 0.001)
        
        if (abs(sigma_g_new - sigma_g) < 0.001 && abs(sigma_e_new - sigma_e) < 0.001) {
            break
        }
        
        sigma_g <- sigma_g_new
        sigma_e <- sigma_e_new
    }
    
    h2 <- sigma_g / (sigma_g + sigma_e)
    cat("Heritability estimate:", h2, "\n")
    
    var_y <- var(y)
    se_h2 <- sqrt((2 * (1 - h2)^2 * h2^2) / n)
    cat("Approx SE:", se_h2, "\n")
    
    write.table(
        data.frame(
            h2 = h2,
            se = se_h2,
            method = "AdjHE",
            n_samples = length(common_ids),
            npc = npc
        ),
        file = out_file,
        row.names = FALSE,
        quote = FALSE
    )
} else {
    stop("Method ", method, " not yet implemented. Use AdjHE.")
}

cat("SNP heritability estimation complete. Results saved to:", out_file, "\n")
