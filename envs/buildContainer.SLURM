#!/bin/bash
#SBATCH --job-name=build       # Name of the job
#SBATCH --output=logs/build.out     # Standard output (%a is the array index)
#SBATCH --error=logs/build.err      # Standard error
#SBATCH --mem=32G                    # Memory per chromosome
#SBATCH --time=03:00:00              # Time limit (HH:MM:SS)

module load apptainer

apptainer build base.sif container.def
