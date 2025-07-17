#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH -p agsmall
#SBATCH -o pcair.out
#SBATCH -e pcair.err
#SBATCH --job-name pcair

WORK=$1           # e.g., /scratch.global/and02709
REF=$2            # unused for now
NAME=$3           # e.g., SMILES_GDA
path_to_repo=$4   # Repo used in other steps not yet performed

module load R/4.4.0-openblas-rocky8
export R_LIBS="/home/gdc/public/Ref/R"

cd $WORK/relatedness
Rscript ${path_to_repo}/src/pca_ir_pipeline.R ${WORK} ${NAME}
