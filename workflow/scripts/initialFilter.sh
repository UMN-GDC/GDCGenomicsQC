#!/bin/bash -l

INPUT=$1
OUTPRE=$2
THREADS=$3
TEMP=$4

mkdir -p $TEMP

INTER_FILEPREFIX=$TEMP/intermediate
INPUT_PREF=$(basename $INPUT)

# Sample missingness initial filter
plink --bfile $INPUT --mind 0.1 --make-bed --out ${INTER_FILEPREFIX}_2  --threads $THREADS

# Marker missingness final filter 
plink --bfile ${INTER_FILEPREFIX}_2 --geno 0.02 --make-bed --out ${INTER_FILEPREFIX}_3  --threads $THREADS

# Sample missingness final filter
plink --bfile ${INTER_FILEPREFIX}_3 --mind 0.02 --make-bed --out $OUTPRE --threads $THREADS

plink --bfile $OUTPRE --indep-pairwise 500 10 0.1 --out $OUTPRE --threads $THREADS
plink --bfile $OUTPRE --extract ${OUTPRE}.prune.in --make-bed --out ${OUTPRE}.LDpruned  --threads $THREADS
