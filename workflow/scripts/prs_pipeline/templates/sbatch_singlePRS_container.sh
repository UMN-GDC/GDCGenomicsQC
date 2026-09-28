#!/bin/bash

prs_pipeline_path=$1
config_path=$2

path_singularity_runner=$(find ${prs_pipeline_path} -name "sandbox_singularity_runner.sh")

sbatch --time=6:00:00 \
  --ntasks=1 \
  --cpus-per-task=16 \
  --mem=64g \
  -o SA_prs_pipeline%j.out \
  --job-name prs_pipeline \
  "${path_singularity_runner}" \
  --C "${config_path}"
