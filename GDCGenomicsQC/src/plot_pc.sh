#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8GB
#SBATCH --time=48:00:00
#SBATCH -p agsmall
#SBATCH -o pc_plot.out
#SBATCH -e pc_plot.err
#SBATCH --job-name pc_plot

WORK=$1
REF=$2
NAME=$3
path_to_repo=$4

module load R

cd $WORK/PCA
Rscript $path_to_repo/src/plot_pca.R $WORK