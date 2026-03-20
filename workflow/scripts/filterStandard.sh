#!/bin/bash -l

INPUT=$1
STAGE=$2
THREADS=$3

mkdir -p $STAGE/intermediates/standard_filter
INTER_FILEPREFIX=$STAGE/intermediates/standard_filter/intermediate

# Minor allele frequency filter
plink2 --pfile $INPUT --freq --out $STAGE/MAF_check --threads $THREADS
plink2 --pfile $INPUT --maf 0.01 --make-pgen --out ${INTER_FILEPREFIX}_6  --threads $THREADS

# Hardy-Weinberg equilibrium check
plink2 --pfile ${INTER_FILEPREFIX}_6 --hardy --out ${INTER_FILEPREFIX}_6  --threads $THREADS
awk '{ if ($9 <0.00001) print $0 }' ${INTER_FILEPREFIX}_6.hwe > ${STAGE}/zoomhwe.hwe
plink2 --pfile ${INTER_FILEPREFIX}_6 --hwe 1e-6 --make-pgen --out ${INTER_FILEPREFIX}_7a  --threads $THREADS
plink2 --pfile ${INTER_FILEPREFIX}_7a --hwe 1e-10 --make-pgen --out ${INTER_FILEPREFIX}_7  --threads $THREADS

# Heterozygosity check:
plink2 --pfile ${INTER_FILEPREFIX}_7 --exclude scripts/inversion.txt --range --indep-pairwise 50 5 0.2 --out $STAGE/indepSNP  --threads $THREADS
plink2 --pfile ${INTER_FILEPREFIX}_7 --extract $STAGE/indepSNP.prune.in --het --out $STAGE/R_check  --threads $THREADS

# This is should be cleaned up
Rscript --no-save scripts/heterozygosity_outliers_list.R
sed 's/"// g' $STAGE/fail-het-qc.txt | awk '{print $1, $2}'> $STAGE/het_fail_ind.txt 
# plink2 --pfile $OUTDIR.QC7 --remove het_fail_ind.txt --make-pgen --out $OUTDIR.QC8
plink2 --pfile ${INTER_FILEPREFIX}_7 --make-pgen --out $STAGE/standardFilter  --threads $THREADS

# Extract snps retained in $OUTDIR.QC8
plink2 --pfile $STAGE/standardFilter --write-snplist --out $STAGE/standardFilter  --threads $THREADS

plink2 --pfile $STAGE/standardFilter --indep-pairwise 500 10 0.1 --out $STAGE/standardFilter  --threads $THREADS
plink2 --pfile $STAGE/standardFilter --extract $STAGE/standardFilter.prune.in --make-pgen --out $STAGE/standardFilter.LDpruned  --threads $THREADS


