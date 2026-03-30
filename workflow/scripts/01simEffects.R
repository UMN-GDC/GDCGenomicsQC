library(argparse)
library(tidyverse)
library(mvtnorm)

parser <- ArgumentParser(description= "Simulate SNP effects controlling for heritability and genetic correlation between ancestries")
parser$add_argument("--ancestry1", type= "character", 
    help = "Filepath to the allele frequencies for ancestry 1 (.afreq file format from Plink)")
parser$add_argument("--ancestry2", type= "character", 
    help = "Filepath to the allele frequencies for ancestry 2 (.afreq file format from Plink)")
parser$add_argument("--seed", type= "integer", default = as.integer(Sys.time()), 
    help = "Specify the desired heritability. Default system time")
parser$add_argument("--heritability", type= "numeric", default = 0.5, 
    help = "Specify the desired heritability. Default 0.5", metavar="h2")
parser$add_argument("--rho", type= "numeric", default = 0.5, metavar="MAF",
    help = "Specify the genetic correlation for each SNP between ancestries. Default 0.5")
parser$add_argument("--maf", type= "numeric", default = 0.05,  metavar= "MAF",
    help = "Specify the minor allele frequency cutoff for causal SNPs. Default 0.05")
parser$add_argument("--selection", type= "numeric", default = 1.0,  metavar= "alpha",
    help = "Specify the type of selection for the causal SNPs. Default 1.0")
args <- parser$parse_args()
#args <- parser$parse_args(c("--ancestry1", "ljkj"))

herit <- args$heritability # might consider having different herits for ancestries
rho <- args$rho   # might consider having different covariance for different ancestries
maf <- args$maf
anc1File <- args$ancestry1
anc2File <- args$ancestry2 # currently only simulates for 2 ancestries can't do a single one right now
selection <- args$selection # currently only simulates for 2 ancestries can't do a single one right now

set.seed(args$seed)

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
combined <- inner_join(anc1, anc2, by = c("#CHROM", "ID", "REF", "ALT"), suffix = c("_EUR", "_AFR"))

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
    matrix(c(combined$sdEff_EUR[i], (rho * herit) / (sqrt(C_EUR) * sqrt(C_AFR)),
      (rho * herit) / (sqrt(C_AFR) * sqrt(C_EUR)), combined$sdEff_AFR[i]), ncol=2, nrow=2))
}

# NOTE: Could induce selection effects
# form the supplemtary note
combined$EUReffs <- betas[,1] * sqrt(2 * combined$ALT_FREQS_EUR * (1- combined$ALT_FREQS_EUR)) ^selection
combined$AFReffs <- betas[,2] * sqrt(2 * combined$ALT_FREQS_AFR * (1- combined$ALT_FREQS_AFR)) ^selection

# Rescale to conserve heritbaility
combined$EUReffs <- combined$EUReffs * sqrt(herit / sum(combined$EUReffs^2))
combined$AFReffs <- combined$AFReffs * sqrt(herit / sum(combined$AFReffs^2))


combined %>%
  select(`#ID` = ID, REF, EUReffs) %>%
  write_delim(paste0(anc1File, "Effects.csv"), delim = "\t")
combined %>%
  select(`#ID` = ID, REF, AFReffs) %>%
  write_delim(paste0(anc2File, "Effects.csv"), delim = "\t")
