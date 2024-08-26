#!/bin/bash

WORK=$1
REF=$2
NAME=$3
path_to_repo=$4
DATATYPE=$5

mkdir $WORK/relatedness
perl $REF/PRIMUS/bin/run_PRIMUS.pl --file ${WORK}/${DATATYPE}/${DATATYPE}.QC8 --genome -t 0.2 -o ${WORK}/relatedness
OUT=$WORK/relatedness/$DATATYPE.QC8_cleaned.genome_maximum_independent_set
# Reformat the unrelated set text file in a suitable format for plink --keep
tail -n +2 "$OUT" > "$OUT.tmp" && mv "$OUT.tmp" "$OUT"
awk '{print "0", $1}' $OUT > "$OUT.tmp" && mv "$OUT.tmp" "$OUT"
# Keep only the unrelated set of individuals determined by PRIMUS
plink --bfile $WORK/$DATATYPE/$DATATYPE.QC8 --keep ${OUT} --output-chr chrMT --make-bed --out $WORK/relatedness/study.$NAME.unrelated
