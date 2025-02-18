#!/bin/bash

WORK=$1
REF=$2
NAME=$3
path_to_repo=$4

mkdir $WORK/PCA
cd $WORK/PCA

#Rscript ${path_to_repo}/src/fraposaRpackage.R

# fraposa operations
# $REF/Fraposa/commvar.sh $REF/PCA_ref/1000G.aligned $WORK/relatedness/study.$NAME.unrelated 1000G.comm study.$NAME.unrelated.comm
$REF/Fraposa/commvar.sh $REF/hgdp1kg.filtered.SNV_INDEL.38.phased.shapeit5/ALL.hgdp1kg.filtered.SNV_INDEL.38.phased.shapeit5 $WORK/relatedness/study.$NAME.unrelated 1000G.comm study.$NAME.unrelated.comm
$REF/Fraposa/fraposa_runner.py --stu_filepref study.$NAME.unrelated.comm 1000G.comm #Main program for Fraposa 
$REF/Fraposa/predstupopu.py 1000G.comm study.$NAME.unrelated.comm 
$REF/Fraposa/plotpcs.py 1000G.comm study.$NAME.unrelated.comm

awk -F '\t' '{print $3}' study.$NAME.unrelated.comm.popu | sort | uniq -c > subpop.txt
awk '{print $1 "\t" $2 "\t" $3}' study.$NAME.unrelated.comm.popu > data.txt
echo "Running subpop.R script"
Rscript ${path_to_repo}/src/subpop.R ${WORK} ${NAME}
rm *.dat
