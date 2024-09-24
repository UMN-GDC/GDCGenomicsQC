#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4GB
#SBATCH --time=2:00:00
#SBATCH -p agsmall
#SBATCH -o rfmix_%a.out
#SBATCH -e rfmix_%a.err
#SBATCH --job-name rfmix

WORK=$1
NAME=$2
CHR=$SLURM_ARRAY_TASK_ID

cd $WORK/rfmix

/home/gdc/and02709/ancestry_OG/rfmix/rfmix \
 -f $WORK/phased/${NAME}.chr${CHR}.phased.vcf.gz \
 -r /home/ljmills/and06841/misc/ALL_phase3_shapeit2_mvncall_integrated_v3plus_nounphased_rsID_genotypes_GRCh38_dbSNP.vcf.gz \
 -m /home/araduski/shared/POPULATION_MAP_FILE.txt \
 -g /home/araduski/shared/genetic_map_hg38.txt \
 -o ancestry_chr${CHR} \
 --chromosome=$CHR
