#!/bin/bash
set -e

# Extract the variants that are commen to both binary PLINK files

REF_PREFIX=$1
SAMP_PREFIX=$2
REF_OUT=$3
SAMP_OUT=$4

cut -f2 ${SAMP_PREFIX}.bim | sort > tmp.rs
plink --bfile ${REF_PREFIX} --extract tmp.rs --out ${REF_OUT} --make-bed
cut -f2 ${REF_OUT}.bim | sort > tmp.rs
plink --bfile ${SAMP_PREFIX} --extract tmp.rs --a1-allele ${REF_OUT}.bim 5 2 --out ${SAMP_OUT} --make-bed
rm tmp.rs

if [ -s "${REF_PREFIX}.popu" ]; then
    cp ${REF_PREFIX}.popu ${REF_OUT}.popu
fi
if [ -s "${SAMP_PREFIX}.popu" ]; then
    cp ${SAMP_PREFIX}.popu ${SAMP_OUT}.popu
fi
