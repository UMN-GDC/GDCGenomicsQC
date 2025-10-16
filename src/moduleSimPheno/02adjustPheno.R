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



# rescale the phenotype to have variance equal to herit
anc1pheno <- read_table(paste0(anc1File, ".sscore"))
anc1SD <- sd(anc1pheno$SCORE1_AVG)
anc1Effects <- read_delim(paste0(anc1File, "Effects.csv"), delim = "\t")
anc1Effects %>%
  mutate(
    EUReffs = EUReffs / anc1SD * sqrt(herit)
  ) %>%
  write_delim(paste0(anc1File, "Effects.csv"), delim = "\t")

anc2pheno <- read_table(paste0(anc2File, ".sscore"))
anc2SD <- sd(anc2pheno$SCORE1_AVG)
anc2Effects <- read_delim(paste0(anc2File, "Effects.csv"), delim = "\t")
anc2Effects %>%
  mutate(
    AFReffs = AFReffs / anc2SD * sqrt(herit)
  ) %>%
  write_delim(paste0(anc2File, "Effects.csv"), delim = "\t")
