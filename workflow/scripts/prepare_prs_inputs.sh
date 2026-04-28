#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  prepare_prs_inputs.sh --sim-dir DIR --out-dir DIR --anc1 AFR --anc2 EUR [--phenotype-index 1] [--gwas-fraction 0.5] [--seed 42] [--plink2-bin plink2]

Creates PRS pipeline inputs from simulated ancestry-specific PLINK files:
  - gwas/target_sumstats.txt and gwas/training_sumstats.txt for PRS-CSx
  - gwas/target_sumstats_singlePRS.txt and gwas/training_sumstats_singlePRS.txt for single-ancestry prs_pipeline
  - anc1_plink_files/<ANC1>_simulation_gwas.{bed,bim,fam}
  - anc1_plink_files/<ANC1>_simulation_study_sample.{bed,bim,fam}
  - anc2_plink_files/<ANC2>_simulation_gwas.{bed,bim,fam}
  - anc2_plink_files/<ANC2>_simulation_study_sample.{bed,bim,fam}
  - metadata/<ANC>_gwas.pheno and metadata/<ANC>_study.pheno
USAGE
}

SIM_DIR=""
OUT_DIR=""
ANC1="AFR"
ANC2="EUR"
PHENO_INDEX="1"
GWAS_FRACTION="0.5"
SEED="42"
PLINK2_BIN="plink2"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sim-dir) SIM_DIR="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --anc1) ANC1="$2"; shift 2 ;;
    --anc2) ANC2="$2"; shift 2 ;;
    --phenotype-index) PHENO_INDEX="$2"; shift 2 ;;
    --gwas-fraction) GWAS_FRACTION="$2"; shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;
    --plink2-bin) PLINK2_BIN="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ -z "$SIM_DIR" || -z "$OUT_DIR" ]]; then
  usage >&2
  exit 2
fi

resolve_plink2_bin() {
  local candidate="$1"

  if command -v "$candidate" >/dev/null 2>&1; then
    command -v "$candidate"
    return 0
  fi

  if [[ -x "$candidate" && ! -d "$candidate" ]]; then
    printf '%s\n' "$candidate"
    return 0
  fi

  if [[ -d "$candidate" && -x "$candidate/plink2" ]]; then
    printf '%s\n' "$candidate/plink2"
    return 0
  fi

  return 1
}

PLINK2_BIN_RESOLVED=$(resolve_plink2_bin "$PLINK2_BIN") || {
  echo "plink2 not found. Tried: $PLINK2_BIN" >&2
  echo "Set prsPipeline.path_plink2 in the config, e.g. /projects/standard/gdc/public/plink2" >&2
  exit 127
}

echo "Using PLINK2: $PLINK2_BIN_RESOLVED"

mkdir -p "$OUT_DIR"/{gwas,anc1_plink_files,anc2_plink_files,metadata,logs,tmp}

make_split_files() {
  local fam="$1"
  local prefix="$2"
  local n
  local n_gwas

  n=$(awk 'NF >= 2 {n++} END {print n+0}' "$fam")
  if [[ "$n" -lt 4 ]]; then
    echo "Need at least 4 samples in $fam to create GWAS/study splits; found $n" >&2
    exit 1
  fi

  n_gwas=$(awk -v n="$n" -v f="$GWAS_FRACTION" 'BEGIN {x=int(n*f); if (x < 2) x=2; if (x > n-2) x=n-2; print x}')

  awk -v seed="$SEED" 'BEGIN {srand(seed)} NF >= 2 {print rand() "\t" $1 "\t" $2}' "$fam" \
    | sort -k1,1g \
    | awk -v n_gwas="$n_gwas" -v gwas="$OUT_DIR/metadata/${prefix}_gwas.keep" -v study="$OUT_DIR/metadata/${prefix}_study.keep" \
        'NR <= n_gwas {print $2, $3 > gwas; next} {print $2, $3 > study}'
}

make_pheno_file() {
  local fam="$1"
  local out="$2"
  local pheno_col=$((5 + PHENO_INDEX))

  awk -v col="$pheno_col" 'BEGIN {print "FID\tIID\tPHENO"} NF >= col && $col != "NA" && $col != "-9" {print $1 "\t" $2 "\t" $col}' "$fam" > "$out"
}

run_glm() {
  local bfile="$1"
  local pheno="$2"
  local out_prefix="$3"

  "$PLINK2_BIN_RESOLVED" \
    --bfile "$bfile" \
    --pheno "$pheno" \
    --pheno-name PHENO \
    --glm hide-covar allow-no-covars \
    --out "$out_prefix"
}

format_sumstats() {
  local glm_prefix="$1"
  local out="$2"
  local glm_file

  glm_file=$(find "$(dirname "$glm_prefix")" -maxdepth 1 -type f -name "$(basename "$glm_prefix")*.glm.linear" | head -n 1)
  if [[ -z "$glm_file" ]]; then
    echo "Could not find PLINK2 .glm.linear output for prefix $glm_prefix" >&2
    exit 1
  fi

  awk '
    BEGIN {OFS="\t"}
    NR == 1 {
      for (i = 1; i <= NF; i++) h[$i] = i
      print "SNP", "CHR", "A1", "A2", "BETA", "SE", "P", "N"
      next
    }
    $h["TEST"] == "ADD" && $h["BETA"] != "NA" && $h["SE"] != "NA" {
      ref = $h["REF"]
      alt = $h["ALT"]
      a1 = $h["A1"]
      a2 = (a1 == alt ? ref : alt)
      p = ("P" in h ? $h["P"] : "NA")
      chr = ("#CHROM" in h ? $h["#CHROM"] : $h["CHROM"])
      print $h["ID"], chr, a1, a2, $h["BETA"], $h["SE"], p, $h["OBS_CT"]
    }
  ' "$glm_file" > "$out"
}

format_single_prs_sumstats() {
  local prscsx_sumstats="$1"
  local out="$2"

  awk '
    BEGIN {OFS="\t"}
    NR == 1 {
      print "rsid", "CHR", "A1", "A2", "beta", "beta_se", "p", "N"
      next
    }
    {print $1, $2, $3, $4, $5, $6, $7, $8}
  ' "$prscsx_sumstats" > "$out"
}

process_ancestry() {
  local anc="$1"
  local subdir="$2"
  local sim_prefix="$SIM_DIR/${anc}_simulation"
  local out_prefix="$OUT_DIR/${subdir}/${anc}_simulation"
  local gwas_prefix="$OUT_DIR/${subdir}/${anc}_simulation_gwas"
  local study_prefix="$OUT_DIR/${subdir}/${anc}_simulation_study_sample"
  local gwas_pheno="$OUT_DIR/metadata/${anc}_gwas.pheno"
  local study_pheno="$OUT_DIR/metadata/${anc}_study.pheno"
  local glm_prefix="$OUT_DIR/gwas/${anc}_glm"

  for ext in bed bim fam; do
    [[ -f "${sim_prefix}.${ext}" ]] || { echo "Missing ${sim_prefix}.${ext}" >&2; exit 1; }
  done

  cp "${sim_prefix}.bed" "${out_prefix}.bed"
  cp "${sim_prefix}.bim" "${out_prefix}.bim"
  cp "${sim_prefix}.fam" "${out_prefix}.fam"

  make_split_files "${out_prefix}.fam" "$anc"

  "$PLINK2_BIN_RESOLVED" --bfile "$out_prefix" --keep "$OUT_DIR/metadata/${anc}_gwas.keep" --make-bed --out "$gwas_prefix"
  "$PLINK2_BIN_RESOLVED" --bfile "$out_prefix" --keep "$OUT_DIR/metadata/${anc}_study.keep" --make-bed --out "$study_prefix"

  make_pheno_file "${gwas_prefix}.fam" "$gwas_pheno"
  make_pheno_file "${study_prefix}.fam" "$study_pheno"
  run_glm "$gwas_prefix" "$gwas_pheno" "$glm_prefix"
}

process_ancestry "$ANC1" "anc1_plink_files"
format_sumstats "$OUT_DIR/gwas/${ANC1}_glm" "$OUT_DIR/gwas/target_sumstats.txt"
format_single_prs_sumstats "$OUT_DIR/gwas/target_sumstats.txt" "$OUT_DIR/gwas/target_sumstats_singlePRS.txt"

process_ancestry "$ANC2" "anc2_plink_files"
format_sumstats "$OUT_DIR/gwas/${ANC2}_glm" "$OUT_DIR/gwas/training_sumstats.txt"
format_single_prs_sumstats "$OUT_DIR/gwas/training_sumstats.txt" "$OUT_DIR/gwas/training_sumstats_singlePRS.txt"

cat > "$OUT_DIR/prs_inputs.env" <<EOF
path_data_root="$OUT_DIR"
target_sumstats_file="$OUT_DIR/gwas/target_sumstats.txt"
training_sumstats_file="$OUT_DIR/gwas/training_sumstats.txt"
target_single_prs_sumstats_file="$OUT_DIR/gwas/target_sumstats_singlePRS.txt"
training_single_prs_sumstats_file="$OUT_DIR/gwas/training_sumstats_singlePRS.txt"
target_gwas_pheno_file="$OUT_DIR/metadata/${ANC1}_gwas.pheno"
target_study_pheno_file="$OUT_DIR/metadata/${ANC1}_study.pheno"
training_gwas_pheno_file="$OUT_DIR/metadata/${ANC2}_gwas.pheno"
training_study_pheno_file="$OUT_DIR/metadata/${ANC2}_study.pheno"
reference_SNPS_bim="$OUT_DIR/anc1_plink_files/${ANC1}_simulation_study_sample"
study_sample_plink="$OUT_DIR/anc1_plink_files/${ANC1}_simulation_study_sample"
study_sample_plink_anc2="$OUT_DIR/anc2_plink_files/${ANC2}_simulation_study_sample"
EOF

N_GWAS_ANC1=$(awk 'END {print NR}' "$OUT_DIR/metadata/${ANC1}_gwas.keep")

cat > "$OUT_DIR/prs_prscsx_generated.conf" <<EOF
path_code="/projects/standard/gdc/public/prs_methods/scripts/PRScsx"
path_data_root="$OUT_DIR"
path_ref_dir="/projects/standard/gdc/public/prs_methods/ref/ref_PRScsx/1kg_ref"
path_plink2="/projects/standard/gdc/public/plink2"

anc1="$ANC1"
anc2="$ANC2"

target_sumstats_file="\${path_data_root}/gwas/target_sumstats.txt"
training_sumstats_file="\${path_data_root}/gwas/training_sumstats.txt"
output_dir="/scratch.global/saonli/prs_pipeline/generated_${ANC1}_${ANC2}"

reference_SNPS_bim="\${path_data_root}/anc1_plink_files/${ANC1}_simulation_study_sample"
study_sample_plink="\${path_data_root}/anc1_plink_files/${ANC1}_simulation_study_sample"
study_sample_plink_anc2="\${path_data_root}/anc2_plink_files/${ANC2}_simulation_study_sample"

prs_pipeline="/projects/standard/gdc/public/prs_methods/scripts/prs_pipeline"
EOF

cat > "$OUT_DIR/prs_single_ancestry_${ANC1}_generated.conf" <<EOF
path_data="$OUT_DIR"
path_repo="/projects/standard/gdc/public/prs_methods/scripts/prs_pipeline"
path_plink2="/projects/standard/gdc/public/plink2"

summary_stats_file="\${path_data}/gwas/target_sumstats_singlePRS.txt"
bim_file_path="\${path_data}/anc1_plink_files/${ANC1}_simulation_study_sample.bim"
study_sample="\${path_data}/anc1_plink_files/${ANC1}_simulation_study_sample"
output_path="\${path_data}/single_ancestry_${ANC1}"

n_total_gwas=$N_GWAS_ANC1
EOF

echo "PRS inputs written under $OUT_DIR"
echo "PRS-CSx config: $OUT_DIR/prs_prscsx_generated.conf"
echo "Single-ancestry config: $OUT_DIR/prs_single_ancestry_${ANC1}_generated.conf"
