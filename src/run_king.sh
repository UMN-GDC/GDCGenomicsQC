#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH -p agsmall
#SBATCH -o king.out
#SBATCH -e king.err
#SBATCH --job-name king

WORK=$1           # e.g., /scratch.global/and02709
REF=$2            # unused for now
NAME=$3           # e.g., SMILES_GDA
path_to_repo=$4   # Repo used in other steps

# Derived paths
ROOT_DIR=$WORK/relatedness
PLINK_FILE=$WORK/Initial_QC/QC4
KING_REPO=/home/gdc/shared/king

mkdir -p "$ROOT_DIR"
cd "$ROOT_DIR" || exit 1

cp ${PLINK_FILE}.bed kin.bed
cp ${PLINK_FILE}.bim kin.bim
cp ${PLINK_FILE}.fam temp.fam

# Check if all FIDs are 0
all_fids_zero=$(awk '{if ($1 != 0) exit 1}' temp.fam && echo "yes" || echo "no")

if [ "$all_fids_zero" == "yes" ]; then
    echo "All FIDs are zero. Using IID as FID for KING and PLINK."
    awk '{print $2, $1}' temp.fam > original_fid_map.txt
    awk '{print $2,$2,$3,$4,$5,$6}' temp.fam > kin.fam
else
    echo "FIDs are not all zero. Using original FID/IID."
    cp temp.fam kin.fam
fi

# Run KING
$KING_REPO -b kin.bed --kinship --prefix kinships

# Run PLINK --genome for IBD estimates (needed for ibdPlot)
module load plink
plink --bfile kin --genome full --out kinships

# Kinship + IBD Plotting
module load R/4.4.2-openblas-rocky8
Rscript $path_to_repo/src/kinship.R $ROOT_DIR kinships

# Subset unrelated/related samples
plink --bfile kin --remove to_exclude.txt --make-bed --out study.$NAME.unrelated
plink --bfile kin --keep to_exclude.txt --make-bed --out study.$NAME.related

# Restore FIDs if modified
if [ "$all_fids_zero" == "yes" ]; then
    echo "Restoring original FIDs in unrelated .fam"
    unrelated_fam="study.$NAME.unrelated.fam"
    awk 'NR==FNR {map[$1]=$2; next} {if ($2 in map) $1=map[$2]; print}' original_fid_map.txt $unrelated_fam > tmp.fam
    mv tmp.fam $unrelated_fam
    echo "Original FIDs restored."
fi
