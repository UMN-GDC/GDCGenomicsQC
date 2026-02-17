#!/bin/bash
#SBATCH --job-name=SM_Master
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4000       # doesn't take much memory for the scheduler script
#SBATCH --time=16:00:00  # This script needs to persist to schedule all downstream scripts
#SBATCH --output=logs/slurm-%j.out
#SBATCH --error=logs/slurm-%j.err

ENV_NAME="snakemake"

if [ -z "$CONDA_EXE" ]; then
    # If CONDA_EXE is empty, Conda isn't loaded. Try to find it.
    # Check common locations if the user didn't specify one.
    if [ -z "$CONDA_BASE_DIR" ]; then
        POSSIBLE_LOCS=(
            "$HOME/miniconda3"
            "$HOME/anaconda3"
            "$HOME/miniforge3"
            "/opt/conda"
        )
        for loc in "${POSSIBLE_LOCS[@]}"; do
            if [ -d "$loc" ]; then
                CONDA_BASE_DIR="$loc"
                break
            fi
        done
    fi
    
    if [ -z "$CONDA_BASE_DIR" ]; then
        echo "Error: Could not find Conda installation. Please set CONDA_BASE_DIR in the script."
        exit 1
    fi
    
    # Source the shell script. This 'activates' the conda command.
    source "${CONDA_BASE_DIR}/etc/profile.d/conda.sh"
else
    # If Conda is already in the environment (e.g. from a login shell), 
    # extract the base path from the executable command.
    CONDA_BASE_DIR=$(dirname $(dirname "$CONDA_EXE"))
    source "${CONDA_BASE_DIR}/etc/profile.d/conda.sh"
fi

echo "Activating environment: $ENV_NAME"
conda activate "$ENV_NAME" || { echo "Failed to activate environment $ENV_NAME"; exit 1; }

echo "Starting Snakemake..."

# snakemake --profile ../profiles/hpc.yaml \
#           --configfile ../config/config.yaml>
snakemake --profile <path/to/profiles/hpc> \
          --configfile <path/to/config.yaml> \
          --use-conda

snakemake --cores=1 --use-conda \
    --configfile </path/to/config.yaml> --directory </path/to/GDCGenomicsQC/workflow> --snakefile </path/to/GDCGenomicsQC/workflow/Snakefile> \
    --executor slurm \
    --report --report-stylesheet /path/to/GDCGenomicsQC/report/stylesheet.css \
    run_train_gnomix

echo "Job finished."

