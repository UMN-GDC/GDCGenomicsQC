#!/bin/bash
# Single-ancestry PRS-CS using PRScsx.py in single-population mode (PRS-CSx with
# one population + one sst file is equivalent to PRS-CS for that ancestry).
# Sources a KEY="VALUE" config file (same contract as run_PRScsx.sh). The
# reformat/score/R2 steps mirror run_PRScsx.sh for consistency.
set -euo pipefail

config_file=""
prs_pipeline="/projects/standard/gdc/public/prs_methods/scripts/prs_pipeline"
target_sumstats_file=""
output_dir=""
reference_SNPS_bim=""
study_sample_plink=""
path_code="/projects/standard/gdc/public/prs_methods/scripts/PRScsx"
path_ref_dir="/projects/standard/gdc/public/prs_methods/ref/ref_PRScsx/1kg_ref"
path_plink2="plink2"
rscript="Rscript"
path_python=""
anc1="AFR"
seed="42"
n_gwas=""

usage() {
  cat <<EOF
Usage: $0 --c CONFIG_FILE

CONFIG_FILE is a KEY="VALUE" shell file (same contract as run_PRScsx.sh):
  target_sumstats_file   aligned target sumstats (prepare_sumstats.R output)
  study_sample_plink     study-sample PLINK prefix (scoring + allele check)
  reference_SNPS_bim     study-sample bim prefix passed to PRScsx.py
  output_dir             method output base (files go to <output>/prs_pipeline/PRScs)
  path_code              directory containing PRScsx.py
  path_ref_dir           PRS-CSx LD reference dir (ldblk_1kg_<pop>/ + snpinfo_mult_1kg_hm3)
  path_plink2            plink2 binary
  rscript                Rscript binary/path for R2 eval (default Rscript)
  path_python            explicit python with numpy+h5py (default: python on PATH)
  anc1                   population code (default AFR)
  seed                   MCMC seed (default 42)
  n_gwas                 GWAS sample size (default: last column of sumstats)
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --c) config_file="$2"; shift 2; break ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 1 ;;
  esac
done

if [[ -n "$config_file" ]]; then
    if [[ ! -f "$config_file" ]]; then
        echo "ERROR: Configuration file not found: $config_file" >&2
        exit 1
    fi
    echo "Loading configuration from: $config_file"
    source "$config_file"
fi

# ---- Validate required vars ----
required_vars=(path_code path_ref_dir path_plink2 anc1 target_sumstats_file \
  output_dir reference_SNPS_bim study_sample_plink)
for v in "${required_vars[@]}"; do
  if [[ -z "${!v:-}" ]]; then
    echo "ERROR: --$v is required (set in $config_file or override defaults)" >&2
    exit 1
  fi
done

if [[ ! -f "${target_sumstats_file}" ]]; then
    echo "ERROR: target sumstats file does not exist: ${target_sumstats_file}" >&2
    exit 1
fi
if [[ ! -d "${path_code}" ]] || [[ ! -f "${path_code}/PRScsx.py" ]]; then
    echo "ERROR: PRS-CSx code not found at ${path_code}/PRScsx.py. Set single_prscs.path_code or run with --download-software." >&2
    exit 1
fi
if [[ ! -d "${path_ref_dir}" ]]; then
    echo "ERROR: PRS-CSx LD reference dir not found: ${path_ref_dir}. Set single_prscs.ld_ref_dir or provision resources." >&2
    exit 1
fi

# ---- Environment ----
# Prefer an explicit python (container or conda with numpy/h5py); fall back to
# 'python' on PATH, then a global conda env (matches run_PRScsx.sh).
if [[ -n "$path_python" ]]; then
    PYTHON_CMD="$path_python"
else
    PYTHON_CMD="python"
    if ! command -v python >/dev/null 2>&1; then
        if [[ -f /projects/standard/gdc/public/envs/load_miniconda3.sh ]]; then
            source /projects/standard/gdc/public/envs/load_miniconda3.sh
        else
            echo "ERROR: 'python' not found in PATH. Run inside the prsv2_latest.sif container or source a conda env." >&2
            exit 1
        fi
    fi
fi
if ! command -v "${path_plink2}" >/dev/null 2>&1 && [[ ! -x "${path_plink2}" ]]; then
    echo "ERROR: '${path_plink2}' not found. Set prsPipeline.path_plink2 to a plink2 binary." >&2
    exit 1
fi

# ---- Derived files ----
final_output_dir=${output_dir}/prs_pipeline/PRScs
mkdir -p "${final_output_dir}"

# Aligned target sumstats layout: SNP CHR BP A1 A2 beta beta_se P n_eff
# (output of prepare_sumstats.R). Reformat to the PRS-CSx sst schema.
target_sst_file_using="${final_output_dir}/target_sumstats_PRScs.txt"
awk 'NR==1 {print "SNP","A1","A2","BETA","SE"} NR>1 {print $1,$4,$5,$6,$7}' \
  "${target_sumstats_file}" > "${target_sst_file_using}"

if [[ -z "$n_gwas" ]]; then
    n_gwas=$(awk 'NR==2 {print $NF}' "${target_sumstats_file}")
fi

echo "Running PRS-CS (single-pop PRS-CSx):"
echo "${PYTHON_CMD} ${path_code}/PRScsx.py --ref_dir=${path_ref_dir} \
  --bim_prefix=${reference_SNPS_bim} \
  --sst_file=${target_sst_file_using} \
  --n_gwas=${n_gwas} \
  --pop=${anc1} \
  --out_dir=${final_output_dir} \
  --out_name=PRScs \
  --seed=${seed}"

"${PYTHON_CMD}" "${path_code}/PRScsx.py" --ref_dir="${path_ref_dir}" \
  --bim_prefix="${reference_SNPS_bim}" \
  --sst_file="${target_sst_file_using}" \
  --n_gwas="${n_gwas}" \
  --pop="${anc1}" \
  --out_dir="${final_output_dir}" \
  --out_name=PRScs \
  --seed="${seed}"

pushd "${final_output_dir}" >/dev/null
  out_prefix="PRScs_${anc1}_pst_eff_a1_b0.5_phiauto"
  combined_file="PRScs_${anc1}_combined_weights.txt"

  echo -e "SNP\tA1\tBETA" > "$combined_file"
  for chr in {1..22}; do
    if [[ -f "${out_prefix}_chr${chr}.txt" ]]; then
      awk -v OFS="\t" '{print $1, $2, $3, $4, $5, $6}' "${out_prefix}_chr${chr}.txt" >> "$combined_file"
    fi
  done

  "${path_plink2}" --bfile "${study_sample_plink}" --score "$combined_file" 2 4 6 header --out "PRScs_${anc1}_score"

  "${rscript}" "${prs_pipeline}/src/PRS_sscore_to_R2.R" "PRScs_${anc1}_score.sscore"
  mv PRS_sscore_R_sqr.txt "${anc1}_PRS_sscore_Rsqr.txt"
  mv adj_PRS_sscore_Rsqr.txt "${anc1}_adj_PRS_sscore_Rsqr.txt"
popd >/dev/null

echo "PRS-CS single-ancestry complete: ${final_output_dir}"