#!/bin/bash

WORK=$1
REF=$2
NAME=$3
path_to_repo=$4
file_to_use=$5


echo "WORK: $WORK"
echo "REF: $REF"
echo "NAME: $NAME"
echo "path_to_repo: $path_to_repo"
echo "file_to_use: $file_to_use"

mkdir $WORK/lifted
for chr in {1..22} X Y; do plink --bfile ${file_to_use} --chr $chr --output-chr chrMT --make-bed --out $WORK/lifted/study.$NAME.lifted.chr${chr};  done
echo "Deleting extra files"
rm prep1.* prep2.* result1.* result2.* result3.* prep.bed updated.snp updated.position updated.chr

# Using genome harmonizer, update strand orientation and flip alleles according to the reference dataset.
# sbatch --wait ${path_to_repo}/src/harmonizer.job ${WORK} ${NAME}
echo "Begin autosomal harmonization"
mkdir -p $WORK/aligned
sbatch --time 12:00:00 --mem 15GB --array 1-22 --wait -N1 ${path_to_repo}/src/harmonizer_individual.job ${WORK} ${NAME} ${REF}
mkdir -p ${WORK}/logs
mkdir -p ${WORK}/logs/errors
mkdir -p ${WORK}/logs/out
mv ${WORK}/*.out ${WORK}/logs/out/
mv ${WORK}/*.err ${WORK}/logs/errors/ 

# Currently reference dataset does not have chrY for alignment, and ChrX has no match with study data
# Hence, we bring the unaligned ChrX and ChrY to the result folder, i.e. skipping alignment
cp $WORK/lifted/study.${NAME}.lifted.chrX.bed $WORK/aligned/study.${NAME}.lifted.chrX.aligned.bed
cp $WORK/lifted/study.${NAME}.lifted.chrX.bim $WORK/aligned/study.${NAME}.lifted.chrX.aligned.bim
cp $WORK/lifted/study.${NAME}.lifted.chrX.fam $WORK/aligned/study.${NAME}.lifted.chrX.aligned.fam
cp $WORK/lifted/study.${NAME}.lifted.chrY.bed $WORK/aligned/study.${NAME}.lifted.chrY.aligned.bed
cp $WORK/lifted/study.${NAME}.lifted.chrY.bim $WORK/aligned/study.${NAME}.lifted.chrY.aligned.bim
cp $WORK/lifted/study.${NAME}.lifted.chrY.fam $WORK/aligned/study.${NAME}.lifted.chrY.aligned.fam


${path_to_repo}/src/genotype_harmonizer_log_reader.sh ${WORK}/aligned 
## Creates genome_harmonizer_full_log.txt inside of the aligned directory

# Merge chromosomes for this step
cd $WORK/aligned
rm mergelist.txt
for chr in {2..22} X Y; do echo study.$NAME.lifted.chr${chr}.aligned >> mergelist.txt; done
plink --bfile study.$NAME.lifted.chr1.aligned --merge-list mergelist.txt --allow-no-sex --make-bed --out study.$NAME.lifted.aligned1
plink --bfile study.$NAME.lifted.aligned1 --split-x 'hg38' no-fail --make-bed --out study.$NAME.lifted.aligned
