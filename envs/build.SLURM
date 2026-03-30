#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24gb
#SBATCH --time=03:00:00
#SBATCH --partition=msismall
#SBATCH --job-name=buildImages
#SBATCH --output=build.out

module load apptainer
export APPTAINER_TMPDIR=/scratch.global/coffm049/apptainer_build_tmp
export SINGULARITY_TMPDIR=/scratch.global/coffm049/apptainer_build_tmp
export APPTAINER_CACHEDIR=/scratch.global/coffm049/apptainer_build_tmp
export SINGULARITY_CACHEDIR=/scratch.global/coffm049/apptainer_build_tmp

#apptainer build --fakeroot ancNreport.sif ancNreport.def # takes ~ 40 mins
apptainer build --fakeroot genomeUtils.sif genomeUtils.def # takes ~ 40 mins
# apptainer build --fakeroot rfmix.sif rfmix.def
# apptainer build --fakeroot mash.sif mash.def

# 4. Optional: Clean up the temporary build files after success
# rm -rf $BUILD_DIR
