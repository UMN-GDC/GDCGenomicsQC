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
parser$add_argument("--heritability", type= "numeric", 
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


