#!/bin/bash -l

source /projects/standard/gdc/public/envs/load_miniconda3-2.sh
conda activate gdcPipeline

FILE=$1
OUTPUT=$2

mkdir -p $OUTPUT

plink --bfile $FILE --missing --out $OUTPUT/initial

# Marker missingness initial filter
plink --bfile $FILE --geno 0.1 --make-bed --out $OUTPUT/QC1

# Sample missingness initial filter
plink --bfile $OUTPUT/QC1 --mind 0.1 --make-bed --out $OUTPUT/QC2

# Marker missingness final filter 
plink --bfile $OUTPUT/QC2 --geno 0.02 --make-bed --out $OUTPUT/QC3

# Sample missingness final filter
plink --bfile $OUTPUT/QC3 --mind 0.02 --make-bed --out $OUTPUT/final

# filtering for linkage disequilibrium
# [ ] May want to have ifelse statment if data is phased or not
# Then could use pairphase for better estimates
plink --bfile $OUTPUT/final --indep-pairwise 500 10 0.1 --out $OUTPUT/final
plink --bfile $OUTPUT/final --extract $OUTPUT/final.prune.in --make-bed --out $OUTPUT/final.LDpruned

