#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8GB
#SBATCH --time=48:00:00
#SBATCH -p agsmall
#SBATCH -o rfmix_all.out
#SBATCH -e rfmix_all.err
#SBATCH --job-name rfmix_all

WORK=$1
NAME=$2
path_to_repo=$3

mkdir $WORK/rfmix
sbatch --time 18:00:00 --mem 64GB --cpus-per-task 30 --array 1-22 --wait -N1 ${path_to_repo}/src/rfmix_individual.sh ${WORK} ${NAME}

module load R/4.4.0-openblas-rocky8
export R_LIBS_USER="/home/gdc/public/Ref/R"

Rscript ${path_to_repo}/src/gai.R ${WORK} ${NAME}

