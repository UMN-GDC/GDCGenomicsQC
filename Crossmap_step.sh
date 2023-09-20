#!/bin/bash
# Adding in the virual environment where CrossMap is installed
source /home/faird/shared/code/external/envs/miniconda3/load_miniconda3.sh
conda activate GDC_pipeline

module load plink
module load python

# Converting to VCF files so that CrossMap can be used
plink --bfile DATA  --recode vcf --out DATA 

# Calling CrossMap
CrossMap.py vcf hg19ToHg38.over.chain.gz DATA.vcf hg38.fa out.hg38.vcf

# Converting back to bim/fam/bed format for the rest of the steps
plink --vcf out.hg38.vcf --make-bed --out hg38_DATA

