#!/bin/bash

REF_PREFIX=$1
SAMP_PREFIX=$2
REF_OUT=$3
SAMP_OUT=$4

cut -f2 ${SAMP_PREFIX}.bim | sort > ${SAMP_OUT}.tmp.rs
plink --bfile ${REF_PREFIX} --extract ${SAMP_OUT}.tmp.rs --out ${REF_OUT} --make-bed
cut -f2 ${REF_OUT}.bim | sort > ${SAMP_OUT}.tmp.rs
plink --bfile ${SAMP_PREFIX} --extract ${SAMP_OUT}.tmp.rs --a1-allele ${REF_OUT}.bim 5 2 --out ${SAMP_OUT} --make-bed

if [ -s "${REF_PREFIX}.popu" ]; then
    cp ${REF_PREFIX}.popu ${REF_OUT}.popu
fi
if [ -s "${SAMP_PREFIX}.popu" ]; then
    cp ${SAMP_PREFIX}.popu ${SAMP_OUT}.popu
fi
