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
#############
# REQUIREMENTS: plink2 and gcta64


######################
# BEGIN MODULE CODE
######################

# Read allele frequencies
anc1 <- read_table(paste0(anc1File, ".afreq"))
# number of causal snps
C_EUR <- sum(anc1$ALT_FREQS > maf)
anc2 <- read_table(paste0(anc2File, ".afreq"))
C_AFR <- sum(anc2$ALT_FREQS > maf)

# causal snps aren't necessarily shared
# NOTE: we might want to run different scenarios where all SNPs aren't necessarily shared
combined <- inner_join(anc1, anc2, by = c("#CHROM", "ID", "REF", "ALT", "PROVISIONAL_REF?"), suffix = c("_EUR", "_AFR"))

#%% simulate effects 
# From "A new method for multiancestry polygenic prediction improves performance across diverse populations"

# compute their variance
combined <- combined %>%
  mutate(sdEff_EUR = 
    case_when(
      is.na(ALT_FREQS_EUR) ~ 0,
      ALT_FREQS_EUR < maf ~ 0,
      .default = sqrt(herit/ C_EUR)
    ), 
  sdEff_AFR =
    case_when(
      is.na(ALT_FREQS_AFR) ~ 0,
      ALT_FREQS_AFR < maf ~ 0,
      .default = sqrt(herit/ C_AFR)
    )
  ) 

# Sample the SNP effects
betas <- matrix(NA, nrow = nrow(combined), ncol=2)
for (i in 1:nrow(combined)) {
    betas[i,1:2] = rmvnorm(n=1, rep(0, 2), 
    matrix(c(combined$sdEff_EUR[i], (rho * herit) / sqrt(C_EUR * C_AFR),
      (rho * herit) / sqrt(C_AFR * C_EUR), combined$sdEff_AFR[i]), ncol=2, nrow=2))
}

# NOTE: Could induce selection effects
combined$EUReffs <- betas[,1] / sqrt(2 * combined$ALT_FREQS_EUR * (1- combined$ALT_FREQS_EUR))
combined$AFReffs <- betas[,1] / sqrt(2 * combined$ALT_FREQS_AFR * (1- combined$ALT_FREQS_AFR))

combined %>%
  select(`#ID` = ID, REF, EUReffs) %>%
  write_delim(paste0(anc1File, "Effects.csv"), delim = "\t")
combined %>%
  select(`#ID` = ID, REF, AFReffs) %>%
  write_delim(paste0(anc2File, "Effects.csv"), delim = "\t")
