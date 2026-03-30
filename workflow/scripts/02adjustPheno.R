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
args <- parser$parse_args()
#args <- parser$parse_args(c("--ancestry1", "ljkj"))

herit <- args$heritability # might consider having different herits for ancestries
rho <- args$rho   # might consider having different covariance for different ancestries
maf <- args$maf
anc1File <- args$ancestry1
anc2File <- args$ancestry2 # currently only simulates for 2 ancestries can't do a single one right now

set.seed(args$seed)



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
