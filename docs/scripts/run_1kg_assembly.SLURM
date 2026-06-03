#!/bin/bash
#SBATCH --job-name=1kg_assembly
#SBATCH --output=logs/1kg_assembly_%j.log
#SBATCH --error=logs/1kg_assembly_%j.err
#SBATCH --time=72:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4

# Setup environment
source /users/4/coffm049/miniconda3/etc/profile.d/conda.sh
conda activate snakemake

cd GDCGenomicsQC/workflow

# Run each step sequentially
gdcgenomicsqc --configfile ../config_reference.yaml kgMeta -j 4
gdcgenomicsqc --configfile ../config_reference.yaml kgData -j 22
gdcgenomicsqc --configfile ../config_reference.yaml kgAssemble -j 8
