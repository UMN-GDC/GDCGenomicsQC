#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=20GB
#SBATCH --time=24:00:00
#SBATCH -p msismall
#SBATCH -o QC.out
#SBATCH -e QC.err
#SBATCH --job-name QC

module load plink
module load R/4.4.0-openblas-rocky8
export R_LIBS_USER="/home/gdc/public/Ref/R"

FILE=$1
OUTDIR=$2
REF=PLACE

mkdir -p $OUTDIR
cd $OUTDIR

# Marker missingness initial filter
plink --bfile $FILE --geno G1 --make-bed --out $OUTDIR.QC1

# Sample missingness initial filter
plink --bfile $OUTDIR.QC1 --mind M1 --make-bed --out $OUTDIR.QC2

# Marker missingness final filter 
plink --bfile $OUTDIR.QC2 --geno G2 --make-bed --out $OUTDIR.QC3

# Sample missingness final filter
plink --bfile $OUTDIR.QC3 --mind M2 --make-bed --out $OUTDIR.QC4

# Gender check:
plink --bfile $OUTDIR.QC4 --check-sex
grep 'PROBLEM' plink.sexcheck| awk '{print$1,$2}'>sex_discrepancy.txt
plink --bfile $OUTDIR.QC4 --remove sex_discrepancy.txt --make-bed --out $OUTDIR.QC5

# Minor allele frequency filter
plink --bfile $OUTDIR.QC5 --freq --out MAF_check
plink --bfile $OUTDIR.QC5 --maf MAF1 --make-bed --out $OUTDIR.QC6

# Hardy-Weinberg equilibrium check
plink --bfile $OUTDIR.QC6 --hardy
awk '{ if ($9 <0.00001) print $0 }' plink.hwe>plinkzoomhwe.hwe
plink --bfile $OUTDIR.QC6 --hwe HWE1 --make-bed --out $OUTDIR.QC7a
plink --bfile $OUTDIR.QC7a --hwe HWE2 --hwe-all --make-bed --out $OUTDIR.QC7

# Heterozygosity check:
plink --bfile $OUTDIR.QC7 --exclude $REF/data/inversion.txt --range --indep-pairwise 50 5 0.2 --out indepSNP
plink --bfile $OUTDIR.QC7 --extract indepSNP.prune.in --het --out R_check
Rscript --no-save $REF/src/heterozygosity_outliers_list.R 
sed 's/"// g' fail-het-qc.txt | awk '{print$1, $2}'> het_fail_ind.txt 
plink --bfile $OUTDIR.QC7 --remove het_fail_ind.txt --make-bed --out $OUTDIR.QC8
