#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --time=2:00:00
#SBATCH -p agsmall
#SBATCH -o GenotypeHarmonizer.out
#SBATCH -e GenotypeHarmonizer.err
#SBATCH --job-name GenotypeHarmonizer

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
for chr in {1..22} X Y; do 
plink --bfile ${file_to_use} --chr $chr --make-bed --out $WORK/lifted/study.$NAME.lifted.chr${chr};  done
## Removed this from the above plink command to align with reference genome # --output-chr chrMT  
echo "Deleting extra files"
rm prep1.* prep2.* result1.* result2.* result3.* prep.bed updated.snp updated.position updated.chr

# Using genome harmonizer, update strand orientation and flip alleles according to the reference dataset.
# sbatch --wait ${path_to_repo}/src/harmonizer.job ${WORK} ${NAME}
echo "Begin autosomal harmonization"
mkdir -p $WORK/aligned
sbatch --time 24:00:00 --mem 15GB --array 1-22 --wait -N1 ${path_to_repo}/src/harmonizer_individual.job ${WORK} ${NAME} ${REF}
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
temp_1=$(ls study.${NAME}.lifted.chr*.aligned.bim)
array=(${temp_1//.bim/})

# Exclude the first element and write to a file
printf "%s\n" "${array[@]:1}" > mergelist.txt
plink --bfile "${array[0]}" --merge-list mergelist.txt --allow-no-sex --make-bed --out study.$NAME.lifted.aligned1
plink --bfile study.$NAME.lifted.aligned1 --split-x 'hg38' no-fail --make-bed --out study.$NAME.lifted.aligned
