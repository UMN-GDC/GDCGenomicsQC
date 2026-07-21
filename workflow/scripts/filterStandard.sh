#!/bin/bash -l

INPUT=$1
STAGE=$2
THREADS=$3

SCRIPTS_DIR="$(cd "$(dirname "$0")" && pwd)"
mkdir -p $STAGE/intermediates/standard_filter
INTER_FILEPREFIX=$STAGE/intermediates/standard_filter/intermediate

# Minor allele frequency filter
plink2 --pfile $INPUT --freq --out $STAGE/MAF_check --threads $THREADS
plink2 --pfile $INPUT --maf 0.01 --make-pgen --out ${INTER_FILEPREFIX}_6 --threads $THREADS

# Hardy-Weinberg equilibrium check
plink2 --pfile ${INTER_FILEPREFIX}_6 --hardy --out ${INTER_FILEPREFIX}_6 --threads $THREADS
awk '$9 < 1e-5' ${INTER_FILEPREFIX}_6.hardy > $STAGE/zoomhwe.hwe
plink2 --pfile ${INTER_FILEPREFIX}_6 --hwe 1e-6 --make-pgen --out ${INTER_FILEPREFIX}_7a --threads $THREADS
plink2 --pfile ${INTER_FILEPREFIX}_7a --hwe 1e-10 --make-pgen --out ${INTER_FILEPREFIX}_7 --threads $THREADS

# Heterozygosity check
INVERSION="$SCRIPTS_DIR/inversion.txt"
if [ -f "$INVERSION" ]; then
    plink2 --pfile ${INTER_FILEPREFIX}_7 --exclude "$INVERSION" --range --indep-pairwise 50 5 0.2 --out $STAGE/indepSNP --threads $THREADS
else
    plink2 --pfile ${INTER_FILEPREFIX}_7 --indep-pairwise 50 5 0.2 --out $STAGE/indepSNP --threads $THREADS
fi
plink2 --pfile ${INTER_FILEPREFIX}_7 --extract $STAGE/indepSNP.prune.in --het --out $STAGE/R_check --threads $THREADS
if [ ! -f $STAGE/R_check.het ]; then
    touch $STAGE/R_check.het
fi

Rscript --no-save $SCRIPTS_DIR/heterozygosity_outliers_list.R $STAGE/R_check.het $STAGE

if [ -f $STAGE/het_fail_ind.txt ]; then
    sed 's/"//g' $STAGE/het_fail_ind.txt | awk '{print $1, $2}' > $STAGE/het_fail_clean.txt
    plink2 --pfile ${INTER_FILEPREFIX}_7 --remove $STAGE/het_fail_clean.txt --make-pgen --out $STAGE/standardFilter --threads $THREADS
else
    plink2 --pfile ${INTER_FILEPREFIX}_7 --make-pgen --out $STAGE/standardFilter --threads $THREADS
fi

# LD-pruned output
plink2 --pfile $STAGE/standardFilter --indep-pairwise 500 10 0.1 --out $STAGE/standardFilter --threads $THREADS
plink2 --pfile $STAGE/standardFilter --extract $STAGE/standardFilter.prune.in --make-pgen --out $STAGE/standardFilter.LDpruned --threads $THREADS
