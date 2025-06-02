#!/bin/bash

WORK=$1           # e.g., /scratch.global/and02709
REF=$2            # unused for now
NAME=$3           # e.g., SMILES_GDA
path_to_repo=$4   # Repo used in other steps not yet performed
DATATYPE=$5       # e.g., full

# Derived paths
ROOT_DIR=$WORK/relatedness
PLINK_FILE=$WORK/$DATATYPE/${DATATYPE}.QC8
KING_REPO=/home/gdc/shared/king

# Create directory
mkdir -p "$ROOT_DIR"
cd "$ROOT_DIR" || exit 1

cp ${PLINK_FILE}.bed kin.bed
cp ${PLINK_FILE}.bim kin.bim
cp ${PLINK_FILE}.fam temp.fam

awk '{print $2,$2,$3,$4,$5,$6}' temp.fam >> kin.fam

$KING_REPO -b kin.bed --kinship --prefix kinships

module load R/4.4.2-openblas-rocky8
Rscript $path_to_repo/src/kinship.R $ROOT_DIR kinships.kin0

module load plink
plink --bfile kin --remove to_exclude.txt --make-bed --out study.$NAME.unrelated
plink --bfile kin --keep to_exclude.txt --make-bed --out study.$NAME.related