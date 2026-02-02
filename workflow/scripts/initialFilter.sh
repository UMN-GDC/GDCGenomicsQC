#!/bin/bash -l

INPUT=$1
STAGE=$2

mkdir -p $STAGE/intermediates/initial_filter

INTER_FILEPREFIX=$STAGE/intermediates/initial_filter/intermediate
INPUT_PREF=$(basename $INPUT)

# Marker missingness initial filter
plink --bfile $INPUT --missing --out $STAGE/initial
plink --bfile $INPUT --out ${INTER_FILEPREFIX}_1 --geno 0.1 --make-bed

# Sample missingness initial filter
plink --bfile ${INTER_FILEPREFIX}_1 --mind 0.1 --make-bed --out ${INTER_FILEPREFIX}_2

# Marker missingness final filter 
plink --bfile ${INTER_FILEPREFIX}_2 --geno 0.02 --make-bed --out ${INTER_FILEPREFIX}_3

# Sample missingness final filter
plink --bfile ${INTER_FILEPREFIX}_3 --mind 0.02 --make-bed --out $STAGE/initialFilter

plink --bfile $STAGE/initialFilter --indep-pairwise 500 10 0.1 --out $STAGE/initialFilter
plink --bfile $STAGE/initialFilter --extract $STAGE/initialFilter.prune.in --make-bed --out $STAGE/initialFilter.LDpruned
