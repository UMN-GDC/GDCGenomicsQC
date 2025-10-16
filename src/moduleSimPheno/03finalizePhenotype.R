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


# Make Noise (remainder of phenotype variance for desired herit)
anc1pheno <- read_table(paste0(anc1File, ".sscore"))
anc1pheno %>%
  mutate(
    error = rnorm(n= dim(anc1pheno)[1],0, sqrt(1-herit)),
    error = error / sd(error) * sqrt(1-herit),
    pheno = SCORE1_AVG + error) %>%
    select(-error) %>%
  write_delim(paste0(anc1File, ".pheno"), delim = "\t")
anc2pheno <- read_table(paste0(anc2File, ".sscore"))
anc2pheno %>%
  mutate(
    error = rnorm(n= dim(anc2pheno)[1],0, sqrt(1-herit)),
    error = error / sd(error) * sqrt(1-herit),
    pheno = SCORE1_AVG + error) %>%
    select(-error) %>%
  write_delim(paste0(anc2File, ".pheno"), delim = "\t")


