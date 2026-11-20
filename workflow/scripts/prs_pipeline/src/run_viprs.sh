#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=30GB
#SBATCH --time=12:00:00
#SBATCH -p msismall
#SBATCH -o Viprs.out
#SBATCH --job-name Viprs

# Prefer viprs from the container (prsv2_latest.sif); fall back to the host
# viprs_env conda environment if viprs_fit is not already available.
###### FUNCTION ######
generate_viprs_sumstats() {
local path_plink2="$1"
local bfile_input="$2"
local out_name="$3"
local covariate_file="$4"
local pheno_file="$5"

if [[ -n "$pheno_file" ]]; then
  local pheno_arg=(--pheno "${pheno_file}")
else
  local pheno_arg=()
fi

for chr in {1..22}; do
  ${path_plink2} \
    --bfile "${bfile_input}" \
    --covar "${covariate_file}" \
    "${pheno_arg[@]}" \
    --chr $chr \
    --glm hide-covar cols=+a1freq,+nobs \
    --out "${out_name}"_${chr}

  awk '
  BEGIN {
    OFS="\t";
    print "#CHROM","POS","ID","REF","ALT1","A1","A1_FREQ","OBS_CT","BETA","SE","T_STAT","P"
  }
  NR>1 && $10=="ADD" {
    print $1,$2,$3,$4,$5,$7,$9,$11,$12,$13,$14,$15
  }
  ' "${out_name}_${chr}.PHENO1.glm.linear" > "${out_name}_${chr}.PHENO1_corrected.glm.linear"
done
}

# Convert pipeline-standard GWAS sumstats (SNP CHR BP A1 A2 beta beta_se P n_eff)
# to the magenpy-standard layout (SNP CHR POS A1 A2 BETA SE PVAL N) by renaming the
# header cells according to the column mapper (e.g. "BP=POS,beta=BETA,..."). Column
# order is preserved and any unrecognized columns are left untouched.
convert_pipeline_sumstats() {
local input_file="$1"
local output_file="$2"
local mapper="$3"
local sep="$4"
awk -F "$sep" -v m="$mapper" '
BEGIN {
    OFS = "\t"
    n = split(m, pairs, ",")
    for (i = 1; i <= n; i++) {
        split(pairs[i], kv, "=")
        rename[kv[1]] = kv[2]
    }
}
NR == 1 {
    for (i = 1; i <= NF; i++) {
        if (i > 1) printf "%s", OFS
        if (($i in rename)) printf "%s", rename[$i]
        else printf "%s", $i
    }
    printf "\n"
    next
}
{
    for (i = 1; i <= NF; i++) {
        if (i > 1) printf "%s", OFS
        printf "%s", $i
    }
    printf "\n"
}
' "$input_file" > "$output_file"
}

usage() {
  cat <<EOF
Usage: $0 --c PATH

Runs VIPRS (GWAS sumstats -> model fit -> score -> evaluate) inside the
prsv2_latest.sif container. Requires a shell-readable config file.

Two mutually-exclusive input modes are supported:

  1. Provided summary statistics (default when bfile_gwas_input is empty):
     provided_sumstats        Path/glob to GWAS sumstats (single file or chr_*).
     viprs_sumstats_format    Format for the sumstats (default: custom).
                              'custom' means the pipeline-standard layout
                              (SNP CHR BP A1 A2 beta beta_se P n_eff); it is
                              converted to the native magenpy layout via
                              viprs_custom_sumstats_mapper. Any native viprs
                              format (plink, plink2, magenpy, fastgwa, ...) is
                              passed through unchanged.
     viprs_custom_sumstats_mapper  Column rename map used to convert the
                              'custom' (pipeline-standard) layout to magenpy,
                              e.g. 'BP=POS,beta=BETA,beta_se=SE,P=PVAL,n_eff=N'.
     viprs_custom_sumstats_sep     Delimiter of the 'custom' input (default: tab).
     viprs_gwas_sample_size   Overall GWAS sample size (used if N not in the file).

  2. Genotype-based sumstats generation (only if bfile_gwas_input is set):
     bfile_gwas_input         PLINK prefix used to generate GWAS sumstats (plink2 --glm).
     covariate_file_gwas      Covariate file for the GWAS (FID IID covariates).

Shared config variables:
  path_data                Root directory for genomic data
  out_path                 Base output directory
  path_plink2              plink2 executable (default: plink2)
  bfile_study_sample       PLINK prefix for the study/sample cohort to score
  covariate_file_study_sample  Covariate file for evaluation (FID IID covariates)
  viprs_ref_glob           Glob for VIPRS LD reference panels (e.g. .../ref_viprs/AFR/chr_*)
  pheno_file               Optional FID/IID/phenotype file passed to plink2 --glm (per-phenotype runs)
  phenotype_file           Optional phenotype file for viprs_evaluate; defaults to study_sample .fam column 6
EOF
}

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

####### Variables #######
path_data=/projects/standard/gdc/public/prs_methods/data/simulated_1000G
out_path=/projects/standard/gdc/public/prs_methods/data/simulated_1000G
path_plink2=plink2
bfile_gwas_input=${path_data}/anc1_plink_files/archived/AFR_simulation_gwas
bfile_study_sample=${path_data}/anc1_plink_files/AFR_simulation_study_sample
covariate_file_gwas=/projects/standard/gdc/public/prs_methods/data/simulated_1000G/prs_pipeline/viprs/gwas/temp/viprs_summary_stats_covar_sex_no_header.txt
covariate_file_study_sample=${path_data}/prs_pipeline/viprs/study_sample_covar.txt
viprs_ref_glob="/projects/standard/gdc/public/prs_methods/ref/ref_viprs/AFR/chr_*"
pheno_file=""              # Optional: FID/IID/phenotype file for plink2 --glm (per-phenotype runs)
phenotype_file=""          # Optional: phenotype file for viprs_evaluate; default = study_sample .fam col 6
viprs_output_dir=""        # Optional: override output directory; default = ${out_path}/prs_pipeline/viprs

## Provided-summary-statistics mode (used when bfile_gwas_input is empty)
provided_sumstats=""          # Path/glob to GWAS sumstats (single file or chr_*)
viprs_sumstats_format="custom"      # viprs_fit --sumstats-format; 'custom' = pipeline-standard layout (SNP CHR BP A1 A2 beta beta_se P n_eff), converted to magenpy before fitting
viprs_custom_sumstats_mapper="BP=POS,beta=BETA,beta_se=SE,P=PVAL,n_eff=N"  # Column rename map used to convert the pipeline-standard layout to magenpy-standard column names
viprs_custom_sumstats_sep="\t"      # Delimiter of the custom (pipeline-standard) sumstats file
viprs_gwas_sample_size=""           # Optional: overall GWAS sample size (--gwas-sample-size)


## If the user doesn't provide a covariate file
# awk 'BEGIN{OFS="\t"; print "FID","IID","SEX"} {print $1,$2,$5}' ${bfile_gwas_input}.fam > ${out_file}_covar_sex.txt

#### # ---- Parse args ----
# Check if the number of arguments is 0
if [ "$#" -eq 0 ]; then
    echo "Error: No config file provided."
    echo "Usage: $0 <config>"
    exit 1
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    --c) config_file="$2"; shift 2; break ;;
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

# ---- Load environment ----
if ! command -v viprs_fit >/dev/null 2>&1; then
    if [[ -f /projects/standard/gdc/public/envs/load_miniconda3.sh ]]; then
        source /projects/standard/gdc/public/envs/load_miniconda3.sh
        conda activate viprs_env
    fi
fi
if ! command -v viprs_fit >/dev/null 2>&1; then
    # Inside prsv2_latest.sif, viprs is installed in the 'multiPRS' micromamba env,
    # which is not on the container's default PATH. Add its bin dir explicitly.
    if [[ -x /opt/conda/envs/multiPRS/bin/viprs_fit ]]; then
        export PATH="/opt/conda/envs/multiPRS/bin:${PATH}"
    fi
fi
if ! command -v viprs_fit >/dev/null 2>&1; then
    echo "ERROR: 'viprs_fit' not found in PATH. Run this inside the prsv2_latest.sif container or activate the viprs_env conda env." >&2
    exit 1
fi

# ---- Output ----
if [[ -n "$viprs_output_dir" ]]; then
  out_dir="${viprs_output_dir}"
else
  out_dir="${out_path}/prs_pipeline/viprs"
fi
mkdir -p "${out_dir}"
out_file=${out_dir}/viprs_summary_stats


# ---- Input dispatch: generate sumstats from genotypes OR use provided sumstats ----
# NOTE: provided_sumstats takes precedence over bfile_gwas_input. This matters when
# running through the pipeline, where bfile_gwas_input may still hold its default value.
if [[ -n "$provided_sumstats" ]]; then
    if [[ -n "$bfile_gwas_input" ]]; then
        echo "WARNING: Both 'bfile_gwas_input' and 'provided_sumstats' are set. Using provided_sumstats." >&2
    fi
    echo "VIPRS input mode: using provided GWAS summary statistics (${provided_sumstats})"
    if [[ "$viprs_sumstats_format" == "custom" ]]; then
        # NOTE: viprs 0.1.3's own 'custom' parser is broken (SumstatsTable.from_file
        # crashes with AttributeError when a parser is supplied without a format name).
        # Instead, convert the pipeline-standard layout to the native 'magenpy' layout
        # via the column mapper, then hand the converted files to --sumstats-format magenpy.
        echo "Converting provided sumstats from pipeline-standard layout to magenpy layout (mapper: ${viprs_custom_sumstats_mapper})"
        converted_dir="${out_dir}/provided_sumstats_converted"
        mkdir -p "${converted_dir}"
        converted_pattern="$(basename "${provided_sumstats}")"
        n_converted=0
        for f in ${provided_sumstats}; do
            [[ -f "$f" ]] || continue
            convert_pipeline_sumstats "$f" "${converted_dir}/$(basename "$f")" "${viprs_custom_sumstats_mapper}" "${viprs_custom_sumstats_sep}"
            n_converted=$((n_converted + 1))
        done
        if [[ $n_converted -eq 0 ]]; then
            echo "ERROR: No provided summary statistics files found matching '${provided_sumstats}'" >&2
            exit 1
        fi
        fit_sumstats="${converted_dir}/${converted_pattern}"
        fit_sumstats_format="magenpy"
    else
        fit_sumstats="${provided_sumstats}"
        fit_sumstats_format="${viprs_sumstats_format}"
    fi
elif [[ -n "$bfile_gwas_input" ]]; then
    echo "VIPRS input mode: generating GWAS summary statistics from genotypes (${bfile_gwas_input})"
    generate_viprs_sumstats ${path_plink2} ${bfile_gwas_input} ${out_file} ${covariate_file_gwas} ${pheno_file}

    mkdir -p ${out_dir}/gwas/logs
    mv ${out_dir}/*.log ${out_dir}/gwas/logs
    mv ${out_dir}/*PHENO* ${out_dir}/gwas
    mkdir -p ${out_dir}/gwas/temp
    mv ${out_dir}/gwas/*PHENO1.* ${out_dir}/gwas/temp

    fit_sumstats="${out_dir}/gwas/viprs_summary_stats_*"
    fit_sumstats_format="plink"
else
    echo "ERROR: No input provided. Set 'bfile_gwas_input' to generate sumstats from genotypes," >&2
    echo "       or set 'provided_sumstats' to a GWAS summary statistics file/glob." >&2
    exit 1
fi

# Build viprs_fit arguments common to both modes
fit_args=(
  -l "${viprs_ref_glob}"
  -s "${fit_sumstats}"
  --sumstats-format "${fit_sumstats_format}"
  --output-dir "${out_dir}"
  --model VIPRS
  --exclude-lrld
  --float-precision float64
)
if [[ -n "$viprs_gwas_sample_size" ]]; then
    fit_args+=(--gwas-sample-size "${viprs_gwas_sample_size}")
fi

viprs_fit "${fit_args[@]}"


echo "Preview of outputed posterior distribution of the variant effect sizes produced from viprs_fit call:"
zcat ${out_dir}/VIPRS_EM.fit.gz | head

viprs_score -f "${out_dir}/VIPRS_EM.fit.gz" \
             --bfile "${bfile_study_sample}" \
             --output-file "${out_dir}/VIPRS_PGS"

echo "Preview of outputed PRS produced from viprs_score call:"
head ${out_dir}/VIPRS_PGS.prs

#### Making the pheno file assuming it is not provided ####
if [[ -n "$phenotype_file" ]]; then
  # Strip header if present so the file matches viprs_evaluate expectations (FID IID phenotype)
  awk 'NR==1 && ($1=="FID" || $1=="fid") {next} {print $1, $2, $3}' OFS="\t" "${phenotype_file}" > "${out_dir}/study_sample_pheno.txt"
else
  awk 'BEGIN{OFS="\t"} {print $1,$2,$6}' ${bfile_study_sample}.fam > "${out_dir}/study_sample_pheno.txt"
fi

# #### Making the testing covariate file if not provided ####
#awk 'BEGIN{OFS="\t"} {print $1,$2,$5}' ${bfile_study_sample}.fam > ${out_dir}/study_sample_covar.txt

viprs_evaluate --prs-file "${out_dir}/VIPRS_PGS.prs" \
                --phenotype-file "${out_dir}/study_sample_pheno.txt" \
                --phenotype-col 2 \
                --covariates-file "${covariate_file_study_sample}" \
                --output-file ${out_dir}/viprs_evaluate_results

echo "Checking results from viprs_evaluate on the test data"
head ${out_dir}/viprs_evaluate_results.eval
