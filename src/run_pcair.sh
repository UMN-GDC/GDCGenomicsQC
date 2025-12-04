#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64GB
#SBATCH --time=18:00:00
#SBATCH -p agsmall
#SBATCH -o pcair.out
#SBATCH -e pcair.err
#SBATCH --job-name pcair

WORK=$1           # e.g., /scratch.global/and02709
REF=$2            # unused for now
NAME=$3           # e.g., SMILES_GDA
path_to_repo=$4   # Repo used in other steps not yet performed

source /projects/standard/gdc/public/envs/load_miniconda3-2.sh
conda activate gdcPipeline

mkdir -p $WORK/pca_ir
cd $WORK/pca_ir
Rscript ${path_to_repo}/src/pca_ir_pipeline.R ${WORK} ${NAME}

Rscript ${path_to_repo}/src/pca_relate_pipeline.R ${WORK} ${NAME}
