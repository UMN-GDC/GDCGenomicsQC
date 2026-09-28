#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24gb
#SBATCH --time=24:00:00
#SBATCH --partition=msismall
#SBATCH --job-name=buildImages
#SBATCH --output=build.out

module load apptainer
mkdir -p /scratch.global/baron063/PRS/apptainer_build_tmp
export APPTAINER_TMPDIR=/scratch.global/baron063/PRS/apptainer_build_tmp
export SINGULARITY_TMPDIR=/scratch.global/baron063/PRS/apptainer_build_tmp
export APPTAINER_CACHEDIR=/scratch.global/baron063/PRS/apptainer_build_tmp
export SINGULARITY_CACHEDIR=/scratch.global/baron063/PRS/apptainer_build_tmp

apptainer build --fakeroot singlePRSHelper.sif singlePRSHelper.def
# apptainer build --fakeroot singlePRS.sif singlePRS.def # takes ~ 40 mins
#apptainer build --fakeroot multiPRS.sif multiPRS.def # takes ~ 180+ mins

# 4. Optional: Clean up the temporary build files after success
# rm -rf $BUILD_DIR
