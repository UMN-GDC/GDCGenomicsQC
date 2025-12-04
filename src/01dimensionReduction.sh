#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --job-name=thinnedAttempt
#SBATCH --ntasks=1
#SBATCH --mem=64g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=coffm049@umn.edu
#SBATCH --output=harmonizationPipeline.out
#SBATCH --error=harmonizationPipeline.err

THREADS=8


#%% UMAP
Rscript 01aUmap.R --eigens testData/1kg.eigenvec --out testData/1kgUmap \
  --npc 50 --neighbors 50 \
  --threads 8 \
  --ncoords 2 \
  --seed $RANDOM

# /projects/standard/gdc/shared/rfmix_ref/ALL_phase3_shapeit2_mvncall_integrated_v3plus_nounphased_rsID_genotypes_GRCh38_dbSNP.vcf.gz

#recode plink bfile to vcf (prepping for popVAE)
plink2 --bfile testData/1kg --recode vcf-iid --out testData/1kg


# Works Run PopVAE using the specific Python from the popvae environment
python /projects/standard/gdc/public/popvae/popvae.py \
    --infile testData/1kg.vcf \
    --out testData/1kgvae \
    --max_epochs 500
