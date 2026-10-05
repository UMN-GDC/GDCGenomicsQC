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
# predlmmAce installs from a private repo, so it needs a token file:
#   apptainer build --fakeroot --build-arg TOKEN_FILE=<path> predlmmAce.sif predlmmAce.def
# see envs/README.md; pass the PATH, never the token itself.
declare -a DEFINITIONS=(
    "ancNreport"
    "genomeUtils"
    "rfmix"
    "mash"
    "phenotypeSim"
    "crossmap"
    "popvae"
    "predlmmAce"
)

# Get the array index
INDEX=$SLURM_ARRAY_TASK_ID
DEF="${DEFINITIONS[$INDEX]}"

echo "Building $DEF..."

BUILD_ARGS=()
if [ "$DEF" = "predlmmAce" ]; then
    # predlmm-ace lives in a private repo; see envs/README.md for token setup.
    TOKEN_FILE="${PREDLMMACE_TOKEN_FILE:-$HOME/.config/gdcgenomicsqc/predlmmace.token}"
    if [ ! -r "$TOKEN_FILE" ]; then
        echo "ERROR: cannot read token file '$TOKEN_FILE'." >&2
        echo "Set PREDLMMACE_TOKEN_FILE or create the default one; see envs/README.md" >&2
        exit 1
    fi
    # Pass the path only: the token itself must never reach --build-arg, since
    # %arguments values are stored in the SIF and recoverable after the push.
    BUILD_ARGS=(--build-arg "TOKEN_FILE=$TOKEN_FILE")
fi

apptainer build --fakeroot "${BUILD_ARGS[@]}" "${DEF}.sif" "${DEF}.def"

echo "Done building $DEF"

# Clean up this task's temp files
rm -rf $APPTAINER_TMPDIR
