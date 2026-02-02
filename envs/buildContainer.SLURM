#!/bin/bash
#SBATCH --job-name=build       
#SBATCH --output=logs/build.out 
#SBATCH --error=logs/build.err   
#SBATCH --mem=32G                 
#SBATCH --time=10:00:00 

module load apptainer

export APPTAINER_CACHEDIR=/scratch.global/coffm049/tmp
export APPTAINER_TMPDIR=/scratch.global/coffm049/tmp
export TMPDIR=/scratch.global/coffm049/tmp

# Create the directories so Apptainer doesn't complain
mkdir -p $APPTAINER_CACHEDIR

apptainer build plink.sif plink.def
apptainer build rfmix.sif rfmix.def

# apptainer build base.sif container.def
