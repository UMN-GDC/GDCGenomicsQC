#!/bin/bash

# Simulate preliminary SNP effects
# ARGS: herit (0-1), SNP effect covariance (0-1), maf (0-1)
#   path to ancestry 1 bim/fam/bed, path 2 ancestry2 bim/fam/bed

ANC1="thinnedAFR"
ANC2="thinnedEUR"
HERIT=0.4
RHO=0.8
MAF=0.05

Rscript 01simEffects.R $HERIT $RHO $MAF $ANC1 $ANC2

# Create phenos first attempt
plink2 --bfile $ANC1 --score ${ANC1}Effects.csv 1 2 3 header --out $ANC1
plink2 --bfile $ANC2 --score ${ANC2}Effects.csv 1 2 3 header --out $ANC2


# Adjust SNP effects
Rscript 02adjustPheno.R $HERIT $RHO $MAF $ANC1 $ANC2

# Generate adjusted phenotype
plink2 --bfile $ANC1 --score ${ANC1}Effects.csv 1 2 3 header --out $ANC1
plink2 --bfile $ANC2 --score ${ANC2}Effects.csv 1 2 3 header --out $ANC2

# finalize phenotype
Rscript 03finalizePhenotype.R $HERIT $RHO $MAF $ANC1 $ANC2



#%% Verify (KODY doesn't need, but might find useful for troubleshooting)

# ESTIMATE SNP heritabilty
# plink2 --bfile $ANC1 --make-grm-bin --out $ANC1
# plink2 --bfile $ANC2 --make-grm-bin --out $ANC2
# # [x] Estimate SNP heritability
# gcta64 --grm $ANC2 --pheno $(ANC2).pheno --reml --mpheno 4 --out $ANC2
# gcta64 --grm $ANC1 --pheno $(ANC1).pheno --reml --mpheno 4 --out $ANC1
# 
# # ESTIMTAE PRS heritability
# plink2 --bfile $ANC1 --pheno $(ANC1).pheno --pheno-name pheno --glm allow-no-covars  --out $ANC1
# plink2 --bfile $ANC2 --pheno $(ANC2).pheno --pheno-name pheno --glm allow-no-covars  --out $ANC2
# 
# plink2 --bfile $ANC2 --score $(ANC2).pheno.glm.linear 3 4 12 header --out $(ANC2)prs
# plink2 --bfile $ANC1 --score $(ANC1).pheno.glm.linear 3 4 12 header --out $(ANC1)prs
