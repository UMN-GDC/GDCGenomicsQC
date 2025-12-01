#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=30GB
#SBATCH --time=1:00:00
#SBATCH -p agsmall
#SBATCH -o QC_initial.out
#SBATCH -e QC_initial.err
#SBATCH --job-name QC_initial

source /home/gdc/public/envs/load_miniconda3-2.sh
conda activate gdcPipeline

FILE=$1
REF=$2

mkdir -p Initial_QC
cd Initial_QC

# Marker missingness initial filter
plink --bfile $FILE --geno 0.1 --make-bed --out QC1

# Sample missingness initial filter
plink --bfile QC1 --mind 0.1 --make-bed --out QC2

# Marker missingness final filter 
plink --bfile QC2 --geno 0.02 --make-bed --out QC3

# Sample missingness final filter
plink --bfile QC3 --mind 0.02 --make-bed --out QC4
plink --bfile QC4 --missing

# filtering for linkage disequilibrium
# [ ] May want to have ifelse statment if data is phased or not
# Then could use pairphase for better estimates
plink --bfile QC4 --indep-pairwise 500 10 0.1 --out QC4
plink --bfile QC4 --extract QC4.prune.in --make-bed --out QC4.LDpruned



