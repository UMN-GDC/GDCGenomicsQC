library(tidyverse)
library(pgenlibr)
library(mvtnorm)

#############
# USER INPUTS

args <- commandArgs(trailingOnly = TRUE)
herit <- as.numeric(args[1]) # might consider having different herits for ancestries
rho <- as.numeric(args[2])   # might consider having different covariance for different ancestries
maf <- as.numeric(args[3])
anc1File <- args[4]
anc2File <- args[5]

eurPRS <- read_table(paste0(anc1File, "prs.sscore"))
eurPheno <- read_table(paste0(anc1File, ".fam"), col_names = c("FID", "IID", "PAT", "MAT", "SEX", "PHENO"))

eurPRS <- eurPRS %>%
  left_join(eurPheno, by = c("#FID" = "FID", "IID"))

# [x] estimate PRS heritability
print(summary(lm(PHENO ~ SCORE1_AVG, eurPRS)))






