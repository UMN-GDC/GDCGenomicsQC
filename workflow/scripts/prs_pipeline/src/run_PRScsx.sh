#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=30GB
#SBATCH --time=5:00:00
#SBATCH -p msismall
#SBATCH -o PRScsx.out
#SBATCH --job-name PRScsx


load_config() {
    local file_path="$1"
    if [[ -f "$file_path" ]]; then
        echo "Loading configuration from: $file_path"
        # Source the file. Variables set in the config will override defaults.
        # Note: 'source' is used for shell-readable KEY="VALUE" files.
        source "$file_path"
        return 0
    else
        echo "ERROR: Configuration file not found at $file_path" >&2
        return 1
    fi
}

usage() {
  cat <<EOF
Usage: $0 \
  --c PATH \
  --path_code PATH \
  --path_data_root PATH \
  --path_ref_dir PATH \
  --path_plink2 PATH \
  --anc1 STR \
  --anc2 STR \
  --target_sumstats_file PATH \
  --training_sumstats_file PATH \
  --output_dir PATH \
  --reference_SNPS_bim PATH \
  --study_sample_plink PATH \
  --study_sample_plink_anc2 PATH

Defaults are shown in the header but script will fail if any are incorrect.
EOF
}

# ---- Defaults (documented only; must be overridden or validated) ----
path_code="/projects/standard/gdc/public/prs_methods/scripts/PRScsx"
path_data_root="/projects/standard/gdc/public/prs_methods/data/simulated_1000G"
path_ref_dir="/projects/standard/gdc/public/prs_methods/ref/ref_PRScsx/1kg_ref"
path_plink2="plink2"
anc1="AFR"
anc2="EUR"
target_sumstats_file="${path_data_root}/gwas/target_sumstats.txt"
training_sumstats_file="${path_data_root}/gwas/training_sumstats.txt"
output_dir="${path_data_root}"
reference_SNPS_bim="${path_data_root}/anc1_plink_files/${anc1}_simulation_study_sample"
study_sample_plink="${path_data_root}/anc1_plink_files/${anc1}_simulation_study_sample"
study_sample_plink_anc2="${path_data_root}/anc2_plink_files/${anc2}_simulation_study_sample"
prs_pipeline="/projects/standard/gdc/public/prs_methods/scripts/prs_pipeline"

# ---- Parse args ----
while [[ $# -gt 0 ]]; do
  case "$1" in
    --c) config_file="$2"; shift 2; break ;;
    --path_code) path_code="$2"; shift 2;;
    --prs_pipeline) prs_pipeline="$2"; shift 2;;
    --path_data_root) path_data_root="$2"; shift 2;;
    --path_ref_dir) path_ref_dir="$2"; shift 2;;
    --path_plink2) path_plink2="$2"; shift 2;;
    --anc1) anc1="$2"; shift 2;;
    --anc2) anc2="$2"; shift 2;;
    --target_sumstats_file) target_sumstats_file="$2"; shift 2;;
    --training_sumstats_file) training_sumstats_file="$2"; shift 2;;
    --output_dir) output_dir="$2"; shift 2;;
    --reference_SNPS_bim) reference_SNPS_bim="$2"; shift 2;;
    --study_sample_plink) study_sample_plink="$2"; shift 2;;
    --study_sample_plink_anc2) study_sample_plink_anc2="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown option: $1"; usage; exit 1;;
  esac
done

if [[ -n "$config_file" ]]; then
    # Case 1: --config was specified. Load it, and it overrides ALL defaults.
    load_config "$config_file" || exit 1
else
    # Case 2: No --config was specified.
    # The variables set by the individual flags in step 3 are used.
    echo "No --c file specified. Using command-line arguments and defaults."
fi

# ---- Validate required args ----
required_vars=(path_code path_data_root path_ref_dir path_plink2 anc1 anc2 \
  target_sumstats_file training_sumstats_file output_dir reference_SNPS_bim \
  study_sample_plink study_sample_plink_anc2)

for v in "${required_vars[@]}"; do
  if [[ -z "${!v:-}" ]]; then
    echo "ERROR: --$v is required"; usage; exit 1
  fi
done

paths_to_check=(target_sumstats_file training_sumstats_file)
for ss_file in ${target_sumstats_file} ${training_sumstats_file}; do
  if [[ ! -f "${ss_file}" ]]; then
    echo "ERROR: --$ss_file does not exist"; usage; exit 1
  fi
done

# ---- Load environment ----
# Prefer python from the container (prsv2_latest.sif); fall back to the host
# conda environment if python is not already available.
if ! command -v python >/dev/null 2>&1; then
    if [[ -f /projects/standard/gdc/public/envs/load_miniconda3.sh ]]; then
        source /projects/standard/gdc/public/envs/load_miniconda3.sh
    else
        echo "ERROR: 'python' not found in PATH. Run this inside the prsv2_latest.sif container or source a conda env." >&2
        exit 1
    fi
fi
if ! command -v "${path_plink2}" >/dev/null 2>&1 && [[ ! -x "${path_plink2}" ]]; then
    echo "ERROR: '${path_plink2}' not found. Set path_plink2 to a plink2 binary (e.g. 'plink2')." >&2
    exit 1
fi

# ---- Derived variables ----
final_output_dir=${output_dir}/prs_pipeline/PRScsx
mkdir -p "${final_output_dir}"

target_sst_file_using="${final_output_dir}/target_sumstats_PRScsx.txt"
training_sst_file_using="${final_output_dir}/training_sumstats_PRScsx.txt"

awk 'BEGIN {OFS="\t"} NR==1 {print "SNP","A1","A2","BETA","SE"} NR>1 {print $1,$3,$4,$5,$6}' \
  "${target_sumstats_file}" > "${target_sst_file_using}"

GWAS_sample_size_target=$(awk 'NR==2 {print $NF}' "${target_sumstats_file}")

awk 'BEGIN {OFS="\t"} NR==1 {print "SNP","A1","A2","BETA","SE"} NR>1 {print $1,$3,$4,$5,$6}' \
  "${training_sumstats_file}" > "${training_sst_file_using}"

GWAS_sample_size_training=$(awk 'NR==2 {print $NF}' "${training_sumstats_file}")

echo "Python call:"
echo "python ${path_code}/PRScsx.py --ref_dir=${path_ref_dir} \
  --bim_prefix=${reference_SNPS_bim} \
  --sst_file=${training_sst_file_using},${target_sst_file_using} \
  --n_gwas=${GWAS_sample_size_training},${GWAS_sample_size_target} \
  --pop=${anc2},${anc1} \
  --out_dir=${final_output_dir} \
  --out_name=PRScsx_joint_${anc1} \
  --seed=42"

python "${path_code}/PRScsx.py" --ref_dir="${path_ref_dir}" \
  --bim_prefix="${reference_SNPS_bim}" \
  --sst_file="${training_sst_file_using},${target_sst_file_using}" \
  --n_gwas="${GWAS_sample_size_training},${GWAS_sample_size_target}" \
  --pop="${anc2},${anc1}" \
  --out_dir="${final_output_dir}" \
  --out_name=PRScsx \
  --seed=42

pushd "${final_output_dir}"
  out_prefix_anc1="PRScsx_${anc1}_pst_eff_a1_b0.5_phiauto"
  combined_file_anc1="PRScsx_${anc1}_combined_weights.txt"

  echo -e "SNP\tA1\tBETA" > "$combined_file_anc1"
  for chr in {1..22}; do
    awk -v OFS="\t" '{print $1, $2, $3, $4, $5, $6}' ${out_prefix_anc1}_chr${chr}.txt >> "$combined_file_anc1"
  done

  out_prefix_anc2="PRScsx_${anc2}_pst_eff_a1_b0.5_phiauto"
  combined_file_anc2="PRScsx_${anc2}_combined_weights.txt"

  echo -e "SNP\tA1\tBETA" > "$combined_file_anc2"
  for chr in {1..22}; do
    awk -v OFS="\t" '{print $1, $2, $3, $4, $5, $6}' ${out_prefix_anc2}_chr${chr}.txt >> "$combined_file_anc2"
  done

  "${path_plink2}" --bfile "${study_sample_plink}" --score "$combined_file_anc1" 2 4 6 header --out PRScsx_joint_${anc1}_score

  Rscript ${prs_pipeline}/src/PRS_sscore_to_R2.R PRScsx_joint_${anc1}_score.sscore
  mv PRS_sscore_R_sqr.txt ${anc1}_PRS_sscore_Rsqr.txt
  mv adj_PRS_sscore_Rsqr.txt ${anc1}_adj_PRS_sscore_Rsqr.txt

  "${path_plink2}" --bfile "${study_sample_plink_anc2}" --score "$combined_file_anc2" 2 4 6 header --out PRScsx_joint_${anc2}_score

  Rscript ${prs_pipeline}/src/PRS_sscore_to_R2.R PRScsx_joint_${anc2}_score.sscore
  mv PRS_sscore_R_sqr.txt ${anc2}_PRS_sscore_Rsqr.txt
  mv adj_PRS_sscore_Rsqr.txt ${anc2}_adj_PRS_sscore_Rsqr.txt
popd

