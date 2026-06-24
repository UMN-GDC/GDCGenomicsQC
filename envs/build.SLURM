#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32gb
#SBATCH --time=01:30:00
#SBATCH --partition=msismall
#SBATCH --job-name=buildImage
#SBATCH --output=build_%a.out
#SBATCH --array=1

# can be array=0-5

module load apptainer

# Unique temp/cache dirs per array task to avoid conflicts
export APPTAINER_TMPDIR=/scratch.global/coffm049/apptainer_build_tmp_${SLURM_ARRAY_TASK_ID}
export SINGULARITY_TMPDIR=/scratch.global/coffm049/apptainer_build_tmp_${SLURM_ARRAY_TASK_ID}
export APPTAINER_CACHEDIR=/scratch.global/coffm049/apptainer_build_cache_${SLURM_ARRAY_TASK_ID}
export SINGULARITY_CACHEDIR=/scratch.global/coffm049/apptainer_build_cache_${SLURM_ARRAY_TASK_ID}

# Create unique directories
mkdir -p $APPTAINER_TMPDIR
mkdir -p $APPTAINER_CACHEDIR

# Array of definitions to build
declare -a DEFINITIONS=(
    "ancNreport"
    "genomeUtils"
    "rfmix"
    "mash"
    "phenotypeSim"
    "crossmap"
)

# Get the array index
INDEX=$SLURM_ARRAY_TASK_ID
DEF="${DEFINITIONS[$INDEX]}"

echo "Building $DEF..."

apptainer build --fakeroot "${DEF}.sif" "${DEF}.def"

echo "Done building $DEF"

# Clean up this task's temp files
rm -rf $APPTAINER_TMPDIR
