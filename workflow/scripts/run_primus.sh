#!/bin/bash

INPUT=$1
OUTPUT=$2
REF=$3
#DATATYPE=$4

plink --bfile ${INPUT} --genome --out ${INPUT}
perl ${REF}/PRIMUS/bin/run_PRIMUS.pl --plink_ibd ${INPUT}.genome -t 0.2 -o $OUTPUT
# perl $REF/PRIMUS/bin/run_PRIMUS.pl --file ${WORK}/${DATATYPE}/${DATATYPE}.QC8 --genome -t 0.2 -o ${WORK}/relatedness # Old technique
# OUT=$WORK/relatedness/$DATATYPE.QC8_cleaned.genome_maximum_independent_set # No longer using their prePRIMUS IBD pipeline!
OUT=$OUTPUT.genome_maximum_independent_set

# Reformat the unrelated set text file in a suitable format for plink --keep
tail -n +2 "$OUT" > "$OUT.tmp" && mv "$OUT.tmp" "$OUT"
awk '{print "0", $1}' $OUT > "$OUT.tmp" && mv "$OUT.tmp" "$OUT"

## Check if the FID and IID are the same in the dataset. If so will need to duplicate the second column into the first column of the outputs of PRIMUS
if diff <(cut -d' ' -f1 $INPUT.fam) <(cut -d' ' -f2 $INPUT.fam) >/dev/null; then
    echo "FID and IID are the same in $INPUT.fam"
    if ! diff <(cut -d' ' -f1 ${OUT}) <(cut -d' ' -f2 ${OUT}) >/dev/null; then
      echo "PRIMUS maximum independent sample has FID and IID does not have the same values so fixing it so that it matches our data using"
      awk '{$1=$2}1' OFS=' ' "${OUT}" > $OUTPUT/temp_file && cp $OUTPUT/temp_file "${OUT}"
    fi
else
    echo "FID and IID are different in $INPUT.fam"
    if diff <(cut -d' ' -f1 ${OUT}) <(cut -d' ' -f2 ${OUT}) >/dev/null; then
      echo "PRIMUS maximum independent sample has FID and IID as the same values so manual fixing is necessary to proceed so that it matches our data using"
      exit 1
    fi
fi

# Keep only the unrelated set of individuals determined by PRIMUS
plink --bfile $INPUT --keep ${OUT} --output-chr chrMT --make-bed --out $OUTPUT/unrelated
