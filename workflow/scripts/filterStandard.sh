#!/bin/bash -l

INPUT=$1
STAGE=$2
THREADS=$3

mkdir -p $STAGE/intermediates/standard_filter
INTER_FILEPREFIX=$STAGE/intermediates/standard_filter/intermediate

# Minor allele frequency filter
plink --bfile $INPUT --freq --out $STAGE/MAF_check --threads $THREADS
plink --bfile $INPUT --maf 0.01 --make-bed --out ${INTER_FILEPREFIX}_6  --threads $THREADS

# Hardy-Weinberg equilibrium check
plink --bfile ${INTER_FILEPREFIX}_6 --hardy --out ${INTER_FILEPREFIX}_6  --threads $THREADS
awk '{ if ($9 <0.00001) print $0 }' ${INTER_FILEPREFIX}_6.hwe > ${STAGE}_zoomhwe.hwe
plink --bfile ${INTER_FILEPREFIX}_6 --hwe 1e-6 --make-bed --out ${INTER_FILEPREFIX}_7a  --threads $THREADS
plink --bfile ${INTER_FILEPREFIX}_7a --hwe 1e-10 --hwe-all --make-bed --out ${INTER_FILEPREFIX}_7  --threads $THREADS

# Heterozygosity check:
plink --bfile ${INTER_FILEPREFIX}_7 --exclude scripts/inversion.txt --range --indep-pairwise 50 5 0.2 --out $STAGE/indepSNP  --threads $THREADS
plink --bfile ${INTER_FILEPREFIX}_7 --extract $STAGE/indepSNP.prune.in --het --out $STAGE/R_check  --threads $THREADS

# This is should be cleaned up
Rscript --no-save scripts/heterozygosity_outliers_list.R
sed 's/"// g' $STAGE/fail-het-qc.txt | awk '{print $1, $2}'> $STAGE/het_fail_ind.txt 
# plink --bfile $OUTDIR.QC7 --remove het_fail_ind.txt --make-bed --out $OUTDIR.QC8
plink --bfile ${INTER_FILEPREFIX}_7 --make-bed --out $STAGE/standardFiltered  --threads $THREADS

# Extract snps retained in $OUTDIR.QC8
plink --bfile $STAGE/standardFiltered --write-snplist --out $STAGE/standardFiltered  --threads $THREADS

plink --bfile $STAGE/standardFiltered --indep-pairwise 500 10 0.1 --out $STAGE/standardFiltered  --threads $THREADS
plink --bfile $STAGE/standardFiltered --extract $STAGE/standardFiltered.prune.in --make-bed --out $STAGE/standardFiltered.LDpruned  --threads $THREADS


