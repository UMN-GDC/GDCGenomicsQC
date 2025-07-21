#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH --time=2:00:00
#SBATCH -p agsmall
#SBATCH -o crossmap.out
#SBATCH -e crossmap.err
#SBATCH --job-name crossmap

WORK=$1
REF=$2
FILE=$3
NAME=$4
path_to_repo=$5


echo "WORK: $WORK"
echo "REF: $REF"
echo "FILE: $FILE"
echo "NAME: $NAME"
echo "path_to_repo: $path_to_repo"

# Since plink denote X chromosome's pseudo-autosomal region as a separate 'XY' chromosome, we want to merge to pass ontto LiftOver/CrossMap. 
# We also reformat the numeric chromsome {1-26} to {1-22, X, Y, MT} for LiftOver/CrossMap
plink --bfile $FILE/$NAME --merge-x no-fail --make-bed --out prep1
plink --bfile prep1 --recode --output-chr 'MT' --out prep2

rm prep.bed updated.snp updated.position updated.chr
awk '{print $1, $4-1, $4, $2}' prep2.map > prep.bed

python ${REF}/CrossMap/CrossMap.py bed ${REF}/CrossMap/GRCh37_to_GRCh38.chain.gz prep.bed study.${NAME}.lifted.bed3

awk '{print $4}' study.$NAME.lifted.bed3 > updated.snp
awk '{print $4, $3}' study.$NAME.lifted.bed3 > updated.position
awk '{print $4, $1}' study.$NAME.lifted.bed3 > updated.chr
plink --file prep2 --extract updated.snp --make-bed --out result1
plink --bfile result1 --update-map updated.position --make-bed --out result2
plink --bfile result2 --update-chr updated.chr --make-bed --out result3
plink --bfile result3 --recode --out study.$NAME.lifted
plink --bfile result3 --recode --make-bed --out study.$NAME.lifted
