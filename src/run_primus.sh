#!/bin/bash

WORK=$1
REF=$2
NAME=$3
path_to_repo=$4
DATATYPE=$5

mkdir -p $WORK/relatedness
plink --bfile ${WORK}/${DATATYPE}/${DATATYPE}.QC8 --genome --out ${WORK}/relatedness/${DATATYPE}.QC8
perl ${REF}/PRIMUS/bin/run_PRIMUS.pl --plink_ibd ${WORK}/relatedness/${DATATYPE}.QC8.genome -t 0.2 -o ${WORK}/relatedness
# perl $REF/PRIMUS/bin/run_PRIMUS.pl --file ${WORK}/${DATATYPE}/${DATATYPE}.QC8 --genome -t 0.2 -o ${WORK}/relatedness # Old technique
# OUT=$WORK/relatedness/$DATATYPE.QC8_cleaned.genome_maximum_independent_set # No longer using their prePRIMUS IBD pipeline!
OUT=$WORK/relatedness/$DATATYPE.QC8.genome_maximum_independent_set

# Reformat the unrelated set text file in a suitable format for plink --keep
tail -n +2 "$OUT" > "$OUT.tmp" && mv "$OUT.tmp" "$OUT"
awk '{print "0", $1}' $OUT > "$OUT.tmp" && mv "$OUT.tmp" "$OUT"

## Check if the FID and IID are the same in the dataset. If so will need to duplicate the second column into the first column of the outputs of PRIMUS
if diff <(cut -d' ' -f1 $WORK/$DATATYPE/$DATATYPE.QC8.fam) <(cut -d' ' -f2 $WORK/$DATATYPE/$DATATYPE.QC8.fam) >/dev/null; then
    echo "FID and IID are the same in ${WORK}/${DATATYPE}/${DATATYPE}.QC8.fam"
    if ! diff <(cut -d' ' -f1 ${OUT}) <(cut -d' ' -f2 ${OUT}) >/dev/null; then
      echo "PRIMUS maximum independent sample has FID and IID does not have the same values so fixing it so that it matches our data using"
      awk '{$1=$2}1' OFS=' ' "${OUT}" > temp_file && cp temp_file "${OUT}"
    fi
else
    echo "FID and IID are different in ${WORK}/${DATATYPE}/${DATATYPE}.QC8.fam"
    if diff <(cut -d' ' -f1 ${OUT}) <(cut -d' ' -f2 ${OUT}) >/dev/null; then
      echo "PRIMUS maximum independent sample has FID and IID as the same values so manual fixing is necessary to proceed so that it matches our data using"
      exit 1
    fi
fi

# Keep only the unrelated set of individuals determined by PRIMUS
plink --bfile $WORK/$DATATYPE/$DATATYPE.QC8 --keep ${OUT} --output-chr chrMT --make-bed --out $WORK/relatedness/study.$NAME.unrelated
