#!/bin/bash -l

INPUT=$1
STAGE=$2
THREADS=$3

mkdir -p $STAGE/intermediates/initial_filter

INTER_FILEPREFIX=$STAGE/intermediates/initial_filter/intermediate
INPUT_PREF=$(basename $INPUT)

# Marker missingness initial filter
plink --bfile $INPUT --missing --out $STAGE/initial --threads $THREADS
plink --bfile $INPUT --out ${INTER_FILEPREFIX}_1 --geno 0.1 --make-bed  --threads $THREADS

# Sample missingness initial filter
plink --bfile ${INTER_FILEPREFIX}_1 --mind 0.1 --make-bed --out ${INTER_FILEPREFIX}_2  --threads $THREADS

# Marker missingness final filter 
plink --bfile ${INTER_FILEPREFIX}_2 --geno 0.02 --make-bed --out ${INTER_FILEPREFIX}_3  --threads $THREADS

# Sample missingness final filter
plink --bfile ${INTER_FILEPREFIX}_3 --mind 0.02 --make-bed --out $STAGE/initialFilter  --threads $THREADS

plink --bfile $STAGE/initialFilter --indep-pairwise 500 10 0.1 --out $STAGE/initialFilter  --threads $THREADS
plink --bfile $STAGE/initialFilter --extract $STAGE/initialFilter.prune.in --make-bed --out $STAGE/initialFilter.LDpruned  --threads $THREADS
