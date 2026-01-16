#!/bin/bash -l

source /projects/standard/gdc/public/envs/load_miniconda3-2.sh
conda activate gdcPipeline

OUT=$1           # e.g., /scratch.global/and02709
INPUT_PREFIX=$2          # e.g., SMILES_GDA
COMB=$3

KING_REPO=/projects/standard/gdc/shared/king

mkdir -p $OUT

BED=${OUT}/kin.bed
BIM=${OUT}/kin.bim
FAM=${OUT}/kin.fam
cp ${INPUT_PREFIX}.bed ${BED}
cp ${INPUT_PREFIX}.bim ${BIM}
cp ${INPUT_PREFIX}.fam ${FAM}

# Check if all FIDs are 0
all_fids_zero=$(awk '{if ($1 != 0) exit 1}' $FAM && echo "yes" || echo "no")

if [ "$all_fids_zero" == "yes" ]; then
    echo "All FIDs are zero. Using IID as FID for KING and PLINK."
    awk '{print $2, $1}' $FAM > $OUT/original_fid_map.txt
    awk '{print $2,$2,$3,$4,$5,$6}' $FAM > $OUT/temp.fam
    mv $OUT/temp.fam $FAM
else
    echo "FIDs are not all zero. Using original FID/IID."
fi

# Run KING
$KING_REPO -b $BED --kinship --prefix $OUT/kinships

# Run PLINK --genome for IBD estimates (needed for ibdPlot)
plink --bfile $OUT/kin --genome $OUT/full --out $OUT/kinships

# Kinship + IBD Plotting
Rscript scripts/kinship.R $OUT $OUT/kinships

# Subset unrelated/related samples
if [ ${COMB} -eq 1 ]; then
    plink --bfile $OUT/kin --make-bed --out $OUT/unrelated
else
    plink --bfile $OUT/kin --remove $OUT/to_exclude.txt --make-bed --out $OUT/unrelated
    plink --bfile $OUT/kin --keep $OUT/to_exclude.txt --make-bed --out $OUT/related
fi

# Restore FIDs if modified
if [ "$all_fids_zero" == "yes" ]; then
    echo "Restoring original FIDs in unrelated .fam"
    unrelated_fam="$OUT/unrelated.fam"
    awk 'NR==FNR {map[$1]=$2; next} {if ($2 in map) $1=map[$2]; print}' $OUT/original_fid_map.txt $unrelated_fam > $OUT/tmp.fam
    mv $OUT/tmp.fam $unrelated_fam
    echo "Original FIDs restored."
fi
