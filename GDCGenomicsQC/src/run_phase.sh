#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8GB
#SBATCH --time=48:00:00
#SBATCH -p agsmall
#SBATCH -o phase_all.out
#SBATCH -e phase_all.err
#SBATCH --job-name phase_all

WORK=$1
REF=$2
NAME=$3
path_to_repo=$4

mkdir $WORK/phased
sbatch --time 4:00:00 --mem 16GB --array 1-22 --wait -N1 ${path_to_repo}/src/phase_individual.sh ${WORK} ${NAME}


