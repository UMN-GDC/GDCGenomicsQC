#!/usr/bin/env bash
#SBATCH --time=12:00:00
#SBATCH --ntasks=16
#SBATCH --mem=100g
#SBATCH --tmp=100g
#SBATCH -p agsmall
#SBATCH -o generate_summary_stats.out
#SBATCH --job-name prs_generate_summary_stats

# set -euo pipefail
# IFS=$'\n\t'

# ========== DEFAULT CONFIGURATION ==========
base_location="/projects/standard/gdc/public/prs_methods/data/simulated_1000G"
#base_location="/projects/standard/gdc/public/prs_methods/data/test/sim_2"
anc1_gwas_input="AFR_simulation_gwas"
anc2_gwas_input="EUR_simulation_gwas"
seed=42
p_pca=/projects/standard/gdc/public/prs_methods/data/simulated_1000G/adjusted_1kgPCs.tsv

path_to_repo="/projects/standard/gdc/public/prs_methods/scripts/prs_pipeline"

### --- FUNCTIONS AND ARGPARSER -----------------------------------------------------------------
usage() {
  cat <<EOF
Usage: $0 [options]

Options:
  -1 <anc1_gwas_input>      Target ancestry plink files gwas split (default: ${anc1_gwas_input})
  -2 <anc2_gwas_input>      Training ancestry plink files gwas split  (default: ${anc2_gwas_input})
  -r <repo_path>            Path to prs_pipeline repo (default: ${path_to_repo})
  -b <base_location>        Full path to workspace (default: ${base_location})
  -S <seed>                 Randomization seed (default:42)
  -p <pca>                 Full path pca file (default: ${p_pca})
  -h                        show this help and exit

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
while getopts ":1:2:r:b:S:p:h" opt; do
  case "$opt" in
    1) anc1_gwas_input="$OPTARG" ;;
    2) anc2_gwas_input="$OPTARG" ;;
    r) path_to_repo="$OPTARG" ;;
    b) base_location="$OPTARG" ;;
    S) seed="$OPTARG" ;;
    p) p_pca="$OPTARG" ;;
    h) usage ;;
    *) usage ;;
  esac
done

# ========== Modules ==========
module load R/4.4.0-openblas-rocky8
module load plink
module load plink/2.00-alpha-091019

REF=/projects/standard/gdc/public/Ref
export R_LIBS_USER=${REF}/R
# ===============================
# FAST GWAS PIPELINE (optimized)
# ===============================

log_file="${base_location}/gen_sum_stats.log"

if [ -f ${log_file} ]; then
  rm ${log_file}
fi

mkdir -p "${base_location}"


# log "STEP 1: Thin dataset for PCA"

# plink --bfile "${base_location}/${anc1_gwas_input}" \
#       --thin ${thin} \
#       --make-bed \
#       --threads 16 \
#       --seed "${seed}" \
#       --out "${base_location}/${anc1_gwas_input}_thin" &>> "${log_file}"

# plink --bfile "${base_location}/${anc2_gwas_input}" \
#       --thin ${thin} \
#       --make-bed \
#       --seed "${seed}" \
#       --threads 16 \
#       --out "${base_location}/${anc2_gwas_input}_thin" &>> "${log_file}"

# ========== PCA ==========
# log "STEP 2: PCA on study sample"
# # USE the plink2 path

# plink2 --bfile "${base_location}/${anc1_gwas_input}_thin" \
#       --pca approx 10 \
#       --threads 2 \
#       --out "${base_location}/${anc1_gwas_input}_pca" &>> "${log_file}"

# plink2 --bfile "${base_location}/${anc2_gwas_input}_thin" \
#       --pca approx 10 \
#       --threads 2 \
#       --out "${base_location}/${anc2_gwas_input}_pca" &>> "${log_file}"

#plink --bfile "${base_location}/${anc1_gwas_input}_thin" \
#      --pca 10 \
#      --threads 2 \
#      --out "${base_location}/${anc1_gwas_input}_pca" &>> "${log_file}"

#plink --bfile "${base_location}/${anc2_gwas_input}_thin" \
#      --pca 10 \
#      --threads 2 \
#      --out "${base_location}/${anc2_gwas_input}_pca" &>> "${log_file}"

# ========== Linear GWAS w/ PCA ==========
log "STEP 1: Running GWAS with PCA covariates"

plink --bfile "${base_location}/${anc1_gwas_input}" \
      --covar "${p_pca}" \
      --linear \
      --allow-no-sex \
      --threads 10 \
      --out "${base_location}/target_sumstats_corrected" \
      &>> "${log_file}"

plink --bfile "${base_location}/${anc2_gwas_input}" \
      --covar "${p_pca}" \
      --linear \
      --allow-no-sex \
      --threads 10 \
      --out "${base_location}/training_sumstats_corrected" \
      &>> "${log_file}"

# ========== Summary stat filtering (1 awk pass) ==========
log "STEP 2: Filtering summary stats"

for type in training target; do
  input="${base_location}/${type}_sumstats_corrected.assoc.linear"
  output="${base_location}/${type}_sumstats_final.txt"

  awk '
    BEGIN {OFS="\t"; print "SNP","CHR","A1","BETA","SE","P","n_eff"}
    NR>1 && $5=="ADD" && $8!=0 && $8!="NA" {
      print $2, $1, $4, $7, $7/$8, $9, $6
    }
  ' "$input" > "$output"
done

# ========== STEP 5: Create target summary stats files ==========
log "STEP 3: Creating target summary stats files via R"

Rscript "${path_to_repo}/src/create_sumstats_files.R" \
  "${base_location}" \
  "target_sumstats_final.txt" \
  "${base_location}/${anc1_gwas_input}.bim" \
  "target_sumstats.txt" \
  "training_sumstats_final.txt" \
  "${base_location}/${anc2_gwas_input}.bim" \
  "training_sumstats.txt" \
  > "${base_location}/create_sumstats_files.log" 2>&1


#Rscript "${path_to_repo}/src/creating_target_sumstats_files_multiple.R" \
#  "${base_location}" \
#  "target_sumstats_filtered_clean.txt" \
#  "${base_location}/${anc1_gwas_input}.bim" \
#  "target_sumstats_prosper.txt" \
#  "training_sumstats_filtered_clean.txt" \
#  "${base_location}/${anc2_gwas_input}.bim" \
#  "training_sumstats_prosper.txt" \
#  > "${base_location}/creating_target_sumstats_files.log" 2>&1


log "Finished generating batch of summary statistic files"
