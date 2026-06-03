#!/bin/bash
#SBATCH --job-name=qc_pipeline
#SBATCH --output=logs/qc_pipeline_%j.log
#SBATCH --error=logs/qc_pipeline_%j.err
#SBATCH --time=72:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4

# Setup environment
source /users/4/coffm049/miniconda3/etc/profile.d/conda.sh
conda activate snakemake

cd GDCGenomicsQC/workflow

# Step 1: Initial QC
gdcgenomicsqc --configfile ../config_qc.yaml full/f1.pgen -j 10

# Step 2: Standard QC
gdcgenomicsqc --configfile ../config_qc.yaml full/f1.b38.f2.pgen -j 10

# Step 3: Ancestry-specific QC
gdcgenomicsqc --configfile ../config_qc.yaml EUR/f1.b38.f2.pgen -j 10
