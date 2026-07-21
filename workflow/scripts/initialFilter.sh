#!/bin/bash -l

INPUT=$1
OUTPRE=$2
THREADS=$3
TEMP=$4

mkdir -p $TEMP

INTER_FILEPREFIX=$TEMP/intermediate
INPUT_PREF=$(basename $INPUT)

# Sample missingness initial filter
plink2 --pfile $INPUT --mind 0.1 --make-pgen --out ${INTER_FILEPREFIX}_2  --threads $THREADS --rm-dup exclude-all \

# Marker missingness final filter 
plink2 --pfile ${INTER_FILEPREFIX}_2 --geno 0.02 --make-pgen --out ${INTER_FILEPREFIX}_3  --threads $THREADS

# Sample missingness final filter
plink2 --pfile ${INTER_FILEPREFIX}_3 --mind 0.02 --make-pgen --out $OUTPRE --threads $THREADS

plink2 --pfile $OUTPRE --indep-pairwise 500 10 0.1 --out $OUTPRE --threads $THREADS
plink2 --pfile $OUTPRE --extract ${OUTPRE}.prune.in --make-pgen --out ${OUTPRE}.LDpruned --threads $THREADS

# Pre-Standard-QC metrics (guaranteed regardless of applyStandardQualityControl)
plink2 --pfile $OUTPRE --freq --out $TEMP/initial_QC --threads $THREADS
plink2 --pfile $OUTPRE --hardy --out $TEMP/initial_QC --threads $THREADS
plink2 --pfile $OUTPRE --indep-pairwise 50 5 0.2 --out $TEMP/het_indep --threads $THREADS
plink2 --pfile $OUTPRE --extract $TEMP/het_indep.prune.in --het --out $TEMP/het_indep --threads $THREADS
[ -f $TEMP/het_indep.het ] || touch $TEMP/het_indep.het
