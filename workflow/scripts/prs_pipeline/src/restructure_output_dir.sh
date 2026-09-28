#!/usr/bin/env bash

# ========== DEFAULT CONFIGURATION ==========
base_location="/projects/standard/gdc/public/prs_methods/data/test/sim_1"
#base_location="/projects/standard/gdc/public/prs_methods/data/test/sim_2"
anc1_prefix="AFR_simulation"
anc2_prefix="EUR_simulation"


### --- FUNCTIONS AND ARGPARSER -----------------------------------------------------------------
usage() {
  cat <<EOF
Usage: $0 [options]

Options:
  -1 <anc1_prefix>      Full path to target ancestry plink files (default: ${anc1_prefix})
  -2 <anc2_prefix>      Full path to training ancestry plink files  (default: ${anc2_prefix})
  -r <repo_path>        Path to prs_pipeline repo (default: ${path_to_repo})
  -b <base_location>    Base location (default: ${base_location})
  -h                    Show this help and exit

Example:
  bash prs_pipeline/src/generate_summary_stat_files.sh -1 AFR_simulation_gwas -2 EUR_simulation_gwas -b /projects/standard/gdc/public/prs_methods/data/test/sim_1

  Using default settings but changing base_location of data
    bash prs_pipeline/src/generate_summary_stat_files.sh -b /projects/standard/gdc/public/prs_methods/data/test/sim_2
EOF
  exit 1
}

log() {
  local msg="$1"
  echo "[$(date '+%F %T')] $msg"
}

### --- PARSE ARGS ------------------------------------------------------------
while getopts ":1:2:r:b:h" opt; do
  case "$opt" in
    1) anc1_prefix="$OPTARG" ;;
    2) anc2_prefix="$OPTARG" ;;
    r) path_to_repo="$OPTARG" ;;
    b) base_location="$OPTARG" ;;
    h) usage ;;
    *) usage ;;
  esac
done

### --- START OF SCRIPT -------------------------------------------------------
mkdir -p "${base_location}/logs"
mv "${base_location}"/*.log "${base_location}/logs"

mkdir -p "${base_location}/gwas"
mv "${base_location}"/*.txt "${base_location}/gwas"
mv "${base_location}"/*sumstats* "${base_location}/gwas"
mv "${base_location}"/ancestry* "${base_location}/gwas"
mv ${base_location}/*pca* ${base_location}/gwas
mv ${base_location}/*assoc* ${base_location}/gwas

mkdir -p "${base_location}/anc1_plink_files/archived"
mv "${base_location}/${anc1_prefix}"* "${base_location}/anc1_plink_files"
mv ${base_location}/randomization_ids_anc1 ${base_location}/anc1_plink_files
mv "${base_location}/anc1_plink_files/"*gwas* "${base_location}/anc1_plink_files/archived"

mkdir -p "${base_location}/anc2_plink_files/archived"
mv "${base_location}/${anc2_prefix}"* "${base_location}/anc2_plink_files"
mv ${base_location}/randomization_ids_anc2 ${base_location}/anc2_plink_files
mv "${base_location}/anc2_plink_files/"*gwas* "${base_location}/anc2_plink_files/archived"
