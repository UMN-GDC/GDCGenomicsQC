library(tidyverse)
eurPRS <- read_table("thinnedAFRprs.sscore")
eurPheno <- read_table("thinnedAFR.fam", col_names = c("FID", "IID", "PAT", "MAT", "SEX", "PHENO"))

eurPRS <- eurPRS %>%
  left_join(eurPheno, by = c("#FID" = "FID", "IID"))

# [x] estimate PRS heritability
summary(lm(PHENO ~ SCORE1_AVG, eurPRS))






