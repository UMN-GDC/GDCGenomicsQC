#!/bin/bash


# conda activate datasci
# sim data
python -m toySim
# Add the hastag to column names
sed -i '1s/^/#/' toySim/toy.vcf

plink2 --vcf toySim/toy.vcf --make-bed --out toySim/toy

rm toySim/toy.vcf toySim/toy.log
