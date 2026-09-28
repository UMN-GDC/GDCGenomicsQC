#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=30GB
#SBATCH --time=12:00:00
#SBATCH -p msismall
#SBATCH -o PROSPER.out
#SBATCH --job-name PROSPER

# Runs PROSPER (multi-ancestry PRS with penalized regression + ensemble learning)
# inside the singleprshelper_latest.sif container. Requires a shell-readable config file.

set -eu

###### FUNCTIONS ######

# Split a PLINK dataset into two subsets (tuning + validation) using random sampling.
# Arguments: $1=bfile_prefix, $2=fraction_for_tuning (0-1), $3=seed, $4=output_dir
split_bfile() {
    local bfile="$1"
    local fraction="$2"
    local seed="$3"
    local out_dir="$4"

    mkdir -p "${out_dir}"

    # Create sample list and shuffle with reproducible seed
    awk '{print $1, $2}' "${bfile}.fam" > "${out_dir}/all_samples.txt"

    # Use awk for reproducible shuffling (avoids dependency on shuf --random-source)
    awk -v s="${seed}" 'BEGIN{srand(s)} {print rand(), $0}' "${out_dir}/all_samples.txt" \
        | sort -n -k1,1 \
        | cut -d' ' -f2- \
        > "${out_dir}/shuffled.txt"

    local n_total n_tuning
    n_total=$(wc -l < "${out_dir}/shuffled.txt")
    n_tuning=$(awk -v n="$n_total" -v f="$fraction" 'BEGIN{printf "%d", n * f + 0.5}')

    if [[ "$n_tuning" -lt 1 ]]; then
        echo "ERROR: Split resulted in 0 tuning samples (n_total=${n_total}, fraction=${fraction})" >&2
        return 1
    fi
    if [[ "$n_tuning" -ge "$n_total" ]]; then
        echo "ERROR: Split resulted in all samples going to tuning (n_total=${n_total}, fraction=${fraction})" >&2
        return 1
    fi

    head -n "$n_tuning" "${out_dir}/shuffled.txt" > "${out_dir}/tuning_keep.txt"
    tail -n +$((n_tuning + 1)) "${out_dir}/shuffled.txt" > "${out_dir}/validation_keep.txt"

    local n_validation
    n_validation=$((n_total - n_tuning))
    echo "  Split ${n_total} samples: ${n_tuning} tuning, ${n_validation} validation"

    # Create PLINK subsets without phenotypes (PROSPER requirement)
    plink2 --bfile "${bfile}" \
        --keep "${out_dir}/tuning_keep.txt" \
        --make-bed \
        --output-missing-phenotype -9 \
        --no-psam-pheno \
        --out "${out_dir}/tuning_bfile" \
        2>/dev/null

    plink2 --bfile "${bfile}" \
        --keep "${out_dir}/validation_keep.txt" \
        --make-bed \
        --output-missing-phenotype -9 \
        --no-psam-pheno \
        --out "${out_dir}/validation_bfile" \
        2>/dev/null

    echo "  PLINK subsets written to ${out_dir}/"
}

usage() {
    cat <<EOF
Usage: $0 --c PATH

Runs PROSPER (lassosum2 -> PROSPER.R -> tuning_testing.R) inside the
singleprshelper_latest.sif container. Requires a shell-readable config file.

Config variables:
  study_sample                 PLINK prefix for the TARGET ancestry study sample
                               (will be split into tuning + validation internally)
  prosper_training_sumstats    GWAS sumstats for the TRAINING/REFERENCE ancestry
                               (format: rsid chr a1 a0 beta beta_se n_eff)
  prosper_target_sumstats      GWAS sumstats for the TARGET ancestry
                               (format: rsid chr a1 a0 beta beta_se n_eff)
                               If empty, uses summary_stats_file from the pipeline.
  prosper_training_anc         Ancestry label for training (e.g. EUR)
  prosper_target_anc           Ancestry label for target (e.g. AFR)
  study_sample_anc2            PLINK prefix for training/reference ancestry samples
  phenotype_info_file          External phenotype file (FID IID phenotype); used instead of .fam column 6
  prosper_tuning_fraction      Fraction of study_sample for tuning (default: 0.7)
  prosper_split_seed           Seed for random sample split (default: 42)
  prosper_ncores               Cores for R scripts (default: 5)
  prosper_package              Path to cloned PROSPER repo
  output_path                  Base output directory
  path_plink2                  plink2 executable (default: plink2)
EOF
}

load_config() {
    local file_path="$1"
    if [[ -f "$file_path" ]]; then
        echo "Loading configuration from: $file_path"
        source "$file_path"
        return 0
    else
        echo "ERROR: Configuration file not found at $file_path" >&2
        return 1
    fi
}

####### DEFAULTS #######
study_sample=""
prosper_training_sumstats=""
prosper_target_sumstats=""
prosper_training_anc=""
prosper_target_anc=""
study_sample_anc2=""
prosper_tuning_fraction=0.7
prosper_split_seed=42
prosper_ncores=5
prosper_package="/projects/standard/gdc/public/prs_methods/scripts/PROSPER"
output_path=""
path_plink2="plink2"
phenotype_info_file=""

# ---- Parse args ----
if [ "$#" -eq 0 ]; then
    echo "Error: No config file provided."
    echo "Usage: $0 --c <config>"
    exit 1
fi

while [[ $# -gt 0 ]]; do
    case "$1" in
        --c) config_file="$2"; shift 2; break ;;
        -h|--help) usage; exit 0;;
        *) echo "Unknown option: $1"; usage; exit 1;;
    esac
done

if [[ -n "${config_file:-}" ]]; then
    load_config "$config_file" || exit 1
else
    echo "ERROR: No --c file specified." >&2
    usage
    exit 1
fi

# ---- Validate required inputs ----
missing=()
[[ -z "$study_sample" ]] && missing+=("study_sample")
[[ -z "$prosper_training_sumstats" ]] && missing+=("prosper_training_sumstats")
[[ -z "$prosper_training_anc" ]] && missing+=("prosper_training_anc")
[[ -z "$prosper_target_anc" ]] && missing+=("prosper_target_anc")
[[ -z "$study_sample_anc2" ]] && missing+=("study_sample_anc2")
if [[ ${#missing[@]} -gt 0 ]]; then
    echo "ERROR: Missing required config variables: ${missing[*]}" >&2
    exit 1
fi

# Determine target sumstats: use prosper_target_sumstats if set, else summary_stats_file
if [[ -z "${prosper_target_sumstats}" ]]; then
    if [[ -n "${summary_stats_file:-}" ]]; then
        prosper_target_sumstats="${summary_stats_file}"
        echo "Using summary_stats_file as target sumstats: ${prosper_target_sumstats}"
    else
        echo "ERROR: Neither prosper_target_sumstats nor summary_stats_file is set." >&2
        exit 1
    fi
fi

# ---- Output directory ----
prosper_dir="${output_path}/prs_pipeline/PROSPER"
mkdir -p "${prosper_dir}/temp"

echo "=========================================="
echo "PROSPER Pipeline"
echo "  Training ancestry: ${prosper_training_anc}"
echo "  Target ancestry:   ${prosper_target_anc}"
echo "  Study sample:      ${study_sample}"
echo "  Training samples:  ${study_sample_anc2}"
echo "  Output:            ${prosper_dir}"
echo "=========================================="

# ---- Step 1: Split study_sample into tuning + validation ----
echo ""
echo "[Step 1] Splitting study_sample into tuning + validation..."
split_bfile "${study_sample}" "${prosper_tuning_fraction}" "${prosper_split_seed}" "${prosper_dir}/temp"

target_tuning_bfile="${prosper_dir}/temp/tuning_bfile"
target_validation_bfile="${prosper_dir}/temp/validation_bfile"

# ---- Step 2: Prepare phenotype files (FID IID phenotype) ----
echo ""
echo "[Step 2] Preparing phenotype files..."

# Helper: extract phenotype for a set of samples from the external phenotype file.
# Falls back to .fam column 6 if the external file has no matches (e.g. training
# ancestry samples not present in target-ancestry phenotype file).
# Usage: extract_pheno <pheno_file> <sample_list_file> <output_file> <fam_file>
extract_pheno() {
    local pheno_file="$1"
    local sample_list="$2"
    local out_file="$3"
    local fam_file="$4"

    if [[ -n "$pheno_file" && -f "$pheno_file" ]]; then
        # Check for header
        local first_field
        first_field=$(head -1 "$pheno_file" | awk '{print $1}')
        local skip=0
        if [[ "$first_field" == "FID" || "$first_field" == "fid" ]]; then
            skip=1
        fi
        awk -v skip="$skip" 'BEGIN{OFS="\t"}
            NR==FNR && NR>skip {pheno[$1,$2]=$3; next}
            ($1,$2) in pheno {print $1, $2, pheno[$1,$2]}
        ' "$pheno_file" "$sample_list" > "$out_file"
    fi

    # If output is empty or pheno_file was missing, fall back to .fam column 6
    if [[ ! -s "$out_file" && -n "$fam_file" && -f "$fam_file" ]]; then
        awk 'BEGIN{OFS="\t"} {print $1, $2, $6}' "$fam_file" > "$out_file"
    fi
}

# Create sample lists (FID IID) for each subset
awk '{print $1, $2}' "${study_sample_anc2}.fam" > "${prosper_dir}/temp/training_samples.txt"

# Training ancestry phenotype (EUR samples — use .fam as fallback since they
# may not appear in the target-ancestry phenotype file)
: > "${prosper_dir}/temp/training_pheno.txt"
extract_pheno "${phenotype_info_file}" \
    "${prosper_dir}/temp/training_samples.txt" \
    "${prosper_dir}/temp/training_pheno.txt" \
    "${study_sample_anc2}.fam"

# Target tuning phenotype (AFR samples — external phenotype file)
: > "${prosper_dir}/temp/target_tuning_pheno.txt"
extract_pheno "${phenotype_info_file}" \
    "${prosper_dir}/temp/tuning_keep.txt" \
    "${prosper_dir}/temp/target_tuning_pheno.txt" \
    "${target_tuning_bfile}.fam"

# Target validation phenotype (AFR samples — external phenotype file)
: > "${prosper_dir}/temp/target_validation_pheno.txt"
extract_pheno "${phenotype_info_file}" \
    "${prosper_dir}/temp/validation_keep.txt" \
    "${prosper_dir}/temp/target_validation_pheno.txt" \
    "${target_validation_bfile}.fam"

echo "  Training pheno:   $(wc -l < "${prosper_dir}/temp/training_pheno.txt") samples"
echo "  Target tuning:    $(wc -l < "${prosper_dir}/temp/target_tuning_pheno.txt") samples"
echo "  Target validation:$(wc -l < "${prosper_dir}/temp/target_validation_pheno.txt") samples"

# ---- Step 3: Convert sumstats to PROSPER format if needed ----
echo ""
echo "[Step 3] Checking sumstats format..."

# PROSPER expects: rsid chr a1 a0 beta beta_se n_eff (tab-delimited, with header)
# Pipeline standard: SNP CHR BP A1 A2 beta beta_se P n_eff

convert_sumstats() {
    local input="$1"
    local output="$2"
    local label="$3"
    local header
    header=$(head -1 "$input")
    if echo "$header" | grep -q "^rsid"; then
        echo "  ${label}: Already in PROSPER format"
        cp "$input" "$output"
    else
        echo "  ${label}: Converting from pipeline-standard to PROSPER format"
        awk 'NR==1 {print "rsid", "chr", "a1", "a0", "beta", "beta_se", "n_eff"; next}
             {print $1, $2, $4, $5, $6, $7, $9}' OFS="\t" "$input" > "$output"
    fi
}

convert_sumstats "${prosper_training_sumstats}" "${prosper_dir}/temp/training_sumstats_prosper.txt" "Training (${prosper_training_anc})"
convert_sumstats "${prosper_target_sumstats}" "${prosper_dir}/temp/target_sumstats_prosper.txt" "Target (${prosper_target_anc})"

# ---- Step 4: Run lassosum2 ----
echo ""
echo "[Step 4] Running lassosum2 for optimal tuning parameters..."

mkdir -p "${prosper_dir}/lassosum2"

Rscript "${prosper_package}/scripts/lassosum2.R" \
    --PATH_package "${prosper_package}" \
    --PATH_out "${prosper_dir}/lassosum2" \
    --PATH_plink "${path_plink2}" \
    --FILE_sst "${prosper_dir}/temp/training_sumstats_prosper.txt,${prosper_dir}/temp/target_sumstats_prosper.txt" \
    --pop "${prosper_training_anc},${prosper_target_anc}" \
    --chrom 1-22 \
    --bfile_tuning "${study_sample_anc2},${target_tuning_bfile}" \
    --pheno_tuning "${prosper_dir}/temp/training_pheno.txt,${prosper_dir}/temp/target_tuning_pheno.txt" \
    --bfile_testing "${study_sample_anc2},${target_validation_bfile}" \
    --pheno_testing "${prosper_dir}/temp/training_pheno.txt,${prosper_dir}/temp/target_validation_pheno.txt" \
    --testing TRUE \
    --cleanup FALSE \
    --NCORES "${prosper_ncores}"

echo "  lassosum2 complete. Optimal params:"
echo "    ${prosper_training_anc}: $(cat "${prosper_dir}/lassosum2/${prosper_training_anc}/optimal_param.txt")"
echo "    ${prosper_target_anc}: $(cat "${prosper_dir}/lassosum2/${prosper_target_anc}/optimal_param.txt")"

# ---- Step 5: Run PROSPER ----
echo ""
echo "[Step 5] Running PROSPER (multi-ancestry penalized regression)..."

mkdir -p "${prosper_dir}/PROSPER"

Rscript "${prosper_package}/scripts/PROSPER.R" \
    --PATH_package "${prosper_package}" \
    --PATH_out "${prosper_dir}/PROSPER" \
    --FILE_sst "${prosper_dir}/temp/training_sumstats_prosper.txt,${prosper_dir}/temp/target_sumstats_prosper.txt" \
    --pop "${prosper_training_anc},${prosper_target_anc}" \
    --lassosum_param "${prosper_dir}/lassosum2/${prosper_training_anc}/optimal_param.txt,${prosper_dir}/lassosum2/${prosper_target_anc}/optimal_param.txt" \
    --chrom 1-22 \
    --NCORES "${prosper_ncores}" \
    --verbose 2

echo "  PROSPER.R complete."

# ---- Step 6: Run tuning_testing (ensemble + evaluation) ----
echo ""
echo "[Step 6] Running tuning_testing (SuperLearner ensemble + evaluation)..."

Rscript "${prosper_package}/scripts/tuning_testing.R" \
    --PATH_plink "${path_plink2}" \
    --PATH_out "${prosper_dir}/PROSPER" \
    --prefix "${prosper_target_anc}" \
    --testing TRUE \
    --bfile_tuning "${target_tuning_bfile}" \
    --pheno_tuning "${prosper_dir}/temp/target_tuning_pheno.txt" \
    --bfile_testing "${target_validation_bfile}" \
    --pheno_testing "${prosper_dir}/temp/target_validation_pheno.txt" \
    --cleanup FALSE \
    --NCORES "${prosper_ncores}" \
    --verbose 2

echo "  tuning_testing complete."

# ---- Summary ----
echo ""
echo "=========================================="
echo "PROSPER Pipeline Complete"
echo "=========================================="
echo "Output directory: ${prosper_dir}/PROSPER/"
echo ""
echo "Key outputs:"
echo "  PRS weights:    ${prosper_dir}/PROSPER/after_ensemble_${prosper_target_anc}/PROSPER_prs_file.txt"
echo "  R2 results:     ${prosper_dir}/PROSPER/after_ensemble_${prosper_target_anc}/R2.txt"
echo "  SL model:       ${prosper_dir}/PROSPER/after_ensemble_${prosper_target_anc}/superlearner_function.RData"
echo ""
if [[ -f "${prosper_dir}/PROSPER/after_ensemble_${prosper_target_anc}/R2.txt" ]]; then
    echo "R2 results:"
    cat "${prosper_dir}/PROSPER/after_ensemble_${prosper_target_anc}/R2.txt"
fi
