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

source /projects/standard/gdc/public/envs/load_miniconda3-2.sh

Rscript ${path_to_repo}/src/gai.R ${WORK} ${NAME}

