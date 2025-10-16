library(tidyverse)
library(pgenlibr)
library(mvtnorm)

#############
# USER INPUTS

herit <- 0.4 # might consider having different herits for ancestries
rho <- 0.8   # might consider having different covariance for different ancestries
maf <- 0.05
#############


# REQUIREMENTS: plink2 and gcta64

################################

system("bash")
# find intersectino of SNPs
awk '{print $2}' ../../data/CTSLEB_sim_data/AFR/AFR_chr1.bim | sort > AFRsnpIDs.txt 
awk '{print $2}' ../../data/CTSLEB_sim_data/EUR/EUR_chr1.bim | sort > EURsnpIDs.txt 


#%% Arbitrarily thin out sample to make it easier for prototyping only using common SNPs
# NOTE: we might want to consider scenarios where not all SNP are shared
plink2 --bfile ../../data/CTSLEB_sim_data/AFR/AFR_chr1 --thin-indiv-count 1000  --thin-count 1000 --make-bed --out thinnedAFR --freq --extract ../../data/CTSLEB_sim_data/EUR/EUR_chr1.bim --maf 0.01
plink2 --bfile ../../data/CTSLEB_sim_data/EUR/EUR_chr1 --thin-indiv-count 1000  --thin-count 1000 --make-bed --out thinnedEUR --freq --extract ../../data/CTSLEB_sim_data/AFR/AFR_chr1.bim --maf 0.01
exit
################################

######################
# BEGIN MODULE CODE
######################
freqFiles <- c("thinnedEUR.afreq", "thinnedAFR.afreq")
.simEffects1 <- function(freqFiles, out, maf = 0.05, rho = 0.8, herit = 0.4) {
  nAnc <- length(freqFiles)
  temp <- freqFiles |>
    # purrr::map(function(x) readr::read_table(x) |> dplyr::mutate(Anc = tools::file_path_sans_ext(x))) |>
    purrr::map(function(x) readr::read_table(x) |> dplyr::mutate(Anc = tools::file_path_sans_ext(x))) |>
    reduce(rbind) |>
    dplyr::mutate(
     C = sum(!is.na(ALT_FREQS)),
     sdEff = dplyr::case_when(
       is.na(ALT_FREQS) ~ 0,
       ALT_FREQS < maf ~ 0,
       .default = sqrt(herit/ C)
     ), .by = Anc 
    ) |>
    select("#CHROM", Anc, ID, sdEff) |>
    tidyr::pivot_wider(id_cols = c("#CHROM", "ID"),
      names_from = Anc, values_from = c("sdEff"))
  # count number of causal snps
  Cs <- colSums(!is.na(temp[,3:ncol(temp)]))

  betas <- matrix(NA, nrow = nrow(temp), ncol=nAnc)
  for (i in 1:nrow(betas)) {
      covar = matrix(NA, nrow = nAnc, ncol = nAnc)
      diag(covar) <- unlist(temp[i, 3:ncol(temp)])

      # replace off diagonal
      for (k in 1:ncol(betas)) {
        for (kprime in 1:ncol(betas)) {
          if (k != kprime) {
              # NOTE: COULD introduce some sort of selection
              # NOTE: COULD introduce difference in rhos 
              entry <- (rho * herit) / sqrt(Cs[k] * Cs[kprime])
              covar[k, kprime] <- entry
              covar[kprime, k] <- entry
        }
      }
    }

      betas[i,] = rmvnorm(n=1, rep(0, nAnc), covar)
  }
  return(betas)
}

betas <- .simEffects1(freqFiles, out="jl", maf = 0.05, rho = 0.8, herit= 0.4)

# [ ] Write generalized saving function
combined %>%
  select(`#ID` = ID, REF, EUReffs) %>%
  write_delim("thinnedEUReffects.csv", delim = "\t")
combined %>%
  select(`#ID` = ID, REF, AFReffs) %>%
  write_delim("thinnedAFReffects.csv", delim = "\t")

# Generate konwn portion of phenotype
system("bash")
plink2 --bfile thinnedAFR --score thinnedAFReffects.csv 1 2 3 header --out thinnedAFR
plink2 --bfile thinnedEUR --score thinnedEUReffects.csv 1 2 3 header --out thinnedEUR
exit

# rescale the phenotype to have variance equal to herit
eurpheno <- read_table("thinnedEUR.sscore")
eurSD <- sd(eurpheno$SCORE1_AVG)
eurEffects <- read_delim("thinnedEUReffects.csv", delim = "\t")
eurEffects %>%
  mutate(
    EUReffs = EUReffs / eurSD * sqrt(herit)
  ) %>%
  write_delim("thinnedEUReffects.csv", delim = "\t")

afrpheno <- read_table("thinnedAFR.sscore")
afrSD <- sd(afrpheno$SCORE1_AVG)
afrEffects <- read_delim("thinnedAFReffects.csv", delim = "\t")
afrEffects %>%
  mutate(
    AFReffs = AFReffs / afrSD * sqrt(herit)
  ) %>%
  write_delim("thinnedAFReffects.csv", delim = "\t")

# Generate konwn portion of phenotype rescaled for controlled heritability
system("bash")
plink2 --bfile thinnedAFR --score thinnedAFReffects.csv 1 2 3 header --out thinnedAFR
plink2 --bfile thinnedEUR --score thinnedEUReffects.csv 1 2 3 header --out thinnedEUR
exit

# Make Noise (remainder of phenotype variance for desired herit)
eurpheno <- read_table("thinnedEUR.sscore")
eurpheno %>%
  mutate(
    error = rnorm(n= dim(eurpheno)[1],0, sqrt(1-herit)),
    error = error / sd(error) * sqrt(1-herit),
    pheno = SCORE1_AVG + error) %>%
    select(-error) %>%
  write_delim("thinnedEUR.pheno", delim = "\t")
afrpheno <- read_table("thinnedAFR.sscore")
afrpheno %>%
  mutate(
    error = rnorm(n= dim(eurpheno)[1],0, sqrt(1-herit)),
    error = error / sd(error) * sqrt(1-herit),
    pheno = SCORE1_AVG + error) %>%
    select(-error) %>%
  write_delim("thinnedAFR.pheno", delim = "\t")


#%% Verify

# [x] make a  GRM
system("bash")
plink2 --bfile thinnedAFR --make-grm-bin --out thinnedAFR
plink2 --bfile thinnedEUR --make-grm-bin --out thinnedEUR
# [x] Estimate SNP heritability
gcta64 --grm thinnedEUR --pheno thinnedEUR.pheno --reml --mpheno 4 --out thinnedEUR
gcta64 --grm thinnedAFR --pheno thinnedAFR.pheno --reml --mpheno 4 --out thinnedAFR

# PRS
plink2 --bfile thinnedAFR --pheno thinnedAFR.pheno --pheno-name pheno --glm allow-no-covars  --out thinnedAFR
plink2 --bfile thinnedEUR --pheno thinnedEUR.pheno --pheno-name pheno --glm allow-no-covars  --out thinnedEUR

plink2 --bfile thinnedEUR --score thinnedEUR.pheno.glm.linear 3 4 12 header --out thinnedEURprs
plink2 --bfile thinnedAFR --score thinnedAFR.pheno.glm.linear 3 4 12 header --out thinnedAFRprs

exit

eurPRS <- read_table("thinnedAFRprs.sscore")
eurPheno <- read_table("thinnedAFR.fam", col_names = c("FID", "IID", "PAT", "MAT", "SEX", "PHENO"))

eurPRS <- eurPRS %>%
  left_join(eurPheno, by = c("#FID" = "FID", "IID"))

# [x] estimate PRS heritability
summary(lm(PHENO ~ SCORE1_AVG, eurPRS))



