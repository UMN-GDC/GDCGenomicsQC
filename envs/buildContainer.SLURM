#!/bin/bash
#SBATCH --job-name=build       # Name of the job
#SBATCH --output=logs/build.out     # Standard output (%a is the array index)
#SBATCH --error=logs/build.err      # Standard error
#SBATCH --mem=32G                    # Memory per chromosome
#SBATCH --time=01:00:00              # Time limit (HH:MM:SS)

module load apptainer

export APPTAINER_CACHEDIR=/scratch.global/coffm049/tmp
export APPTAINER_TMPDIR=/scratch.global/coffm049/tmp

# Create the directories so Apptainer doesn't complain
mkdir -p $APPTAINER_CACHEDIR

apptainer build plink.sif plink.def

# apptainer build base.sif container.def
