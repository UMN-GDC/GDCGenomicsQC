#!/bin/bash 
#SBATCH --time=1:00:00
#SBATCH --ntasks=8
#SBATCH --mem=16g
#SBATCH --tmp=32g
#SBATCH -p agsmall
#SBATCH -o split_plink_data.out
#SBATCH --job-name split_plink_data

# set -euo pipefail
# shopt -s extglob

### --- CONFIGURABLE DEFAULTS -------------------------------------------------
# Defaults (can be overridden via CLI)
plink_file_anc1=""
plink_file_anc2=""

train_percent=50
valid_percent=0
test_percent=50
rand_seed=42
no_plink=0
path_to_repo=/projects/standard/gdc/public/prs_methods/scripts/prs_pipeline # place the code is located
# path_R_packages=/projects/standard/gdc/public/Ref/R
# export R_LIBS_USER="$path_R_packages"


### --- FUNCTIONS ------------------------------------------------------------
usage() {
  cat <<EOF
Usage: $0 [options]

Split PLINK files into train/val/test subsets. When validation percent is 0
(default), produces a 2-way train/test split instead.

Accepts one (-1) or two (-1 and -2) PLINK prefixes. With one prefix, only
that dataset is split. With two, both are split independently (e.g. separate
ancestry batches).

Options:
  -1 <path>    PLINK prefix (target ancestry). Required.
  -2 <path>    Second PLINK prefix (training ancestry). Optional.
  -r <path>    Path to prs_pipeline repo (default: ${path_to_repo})
  -t <pct>     Training percent (default: 50)
  -v <pct>     Validation percent (default: 0; set >0 for 3-way split)
  -T <pct>     Test percent (default: 50)
  -S <int>     Random seed (default: 42)
  -N           Sample lists only (skip PLINK file generation)
  -h           Show this help and exit

Examples:
  Single ancestry (1 input):
    bash $0 -1 /path/to/study_sample

  Two ancestries (2 inputs):
    bash $0 -1 /path/to/AFR -2 /path/to/EUR

  3-way split:
    bash $0 -1 /path/to/study_sample -t 70 -v 15 -T 15

  Sample lists only:
    bash $0 -1 /path/to/study_sample -N
EOF
  exit 1
}

log() {
  local msg="$1"
  echo "[$(date '+%F %T')] $msg"
}



### --- PARSE ARGS ------------------------------------------------------------
while getopts ":1:2:r:t:T:v:S:Nh" opt; do
  case "$opt" in
    1) plink_file_anc1="$OPTARG" ;;
    2) plink_file_anc2="$OPTARG" ;;
    r) path_to_repo="$OPTARG" ;;
    t) train_percent="$OPTARG" ;;
    T) test_percent="$OPTARG" ;;
    v) valid_percent="$OPTARG" ;;
    S) rand_seed="$OPTARG" ;;
    N) no_plink=1 ;;
    h) usage ;;
    *) usage ;;
  esac
done

if [[ -z "$plink_file_anc1" ]]; then
    echo "ERROR: -1 is required"
    usage
fi

# Load environments / modules
# Inside the Singularity container, plink/plink2/python are already in PATH
if [[ -z "${SINGULARITY_CONTAINER:-}" ]]; then
    module load plink
    module load plink/2.00-alpha-091019
fi

base_location=$(dirname "$plink_file_anc1")
mkdir -p "${base_location}/randomization_ids_anc1"
if [[ -n "$plink_file_anc2" ]]; then
    mkdir -p "${base_location}/randomization_ids_anc2"
fi

# Split phenotype / plink files
if [[ -z "${SINGULARITY_CONTAINER:-}" ]]; then
    source /projects/standard/gdc/public/envs/load_miniconda3.sh
fi
log "Splitting plink data for $plink_file_anc1"
run_split() {
    local prefix="$1"
    local out_dir="$2"
    if [ "${no_plink}" -eq 1 ]; then
        python "${path_to_repo}/src/split_plink_samples.py" "${prefix}" \
            --train "${train_percent}" \
            --val "${valid_percent}" \
            --test "${test_percent}" \
            --seed ${rand_seed} \
            --no_plink
    else
        python "${path_to_repo}/src/split_plink_samples.py" "${prefix}" \
            --train "${train_percent}" \
            --val "${valid_percent}" \
            --test "${test_percent}" \
            --seed ${rand_seed}
    fi
    mv "${base_location}"/*samples.txt "${base_location}/${out_dir}"
}

run_split "$plink_file_anc1" randomization_ids_anc1
if [[ -n "$plink_file_anc2" ]]; then
    log "Splitting plink data for $plink_file_anc2"
    run_split "$plink_file_anc2" randomization_ids_anc2
fi

if [[ -z "${SINGULARITY_CONTAINER:-}" ]]; then
    conda deactivate
    conda deactivate
fi


if [[ "${valid_percent}" -gt 0 ]]; then
  echo "Finished generating inputs (3-way split): train=${train_percent}%, val=${valid_percent}%, test=${test_percent}%"
else
  echo "Finished generating inputs (2-way split): train=${train_percent}%, test=${test_percent}%"
fi

