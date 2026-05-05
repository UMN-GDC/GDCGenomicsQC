#!/usr/bin/env bash
set -euo pipefail

METHOD=""
PRS_INPUTS_ENV=""
RESOURCE_DIR=""
OUT_DIR=""
DONE=""
LD_REF_DIR=""
LD_REF_PREFIX=""
LD_MATRIX_DIR=""
SOFTWARE_DIR=""

usage() {
  cat <<'USAGE'
Usage:
  run_prs_pipeline_adapter.sh --method METHOD --prs-inputs-env FILE --resource-dir DIR --out-dir DIR --done FILE [method options]

The adapter validates the GDC-generated PRS inputs and writes a method manifest.
If PRS_METHOD_COMMAND is set, the command is evaluated with common PRS input
variables exported in the environment. If PRS_METHOD_COMMAND is empty, the rule
stages and validates inputs only and writes a READY_NO_COMMAND status.
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --method)
      [[ $# -ge 2 ]] || { echo "Missing value for $1" >&2; exit 2; }
      METHOD="$2"; shift 2 ;;
    --prs-inputs-env)
      [[ $# -ge 2 ]] || { echo "Missing value for $1" >&2; exit 2; }
      PRS_INPUTS_ENV="$2"; shift 2 ;;
    --resource-dir)
      [[ $# -ge 2 ]] || { echo "Missing value for $1" >&2; exit 2; }
      RESOURCE_DIR="$2"; shift 2 ;;
    --out-dir)
      [[ $# -ge 2 ]] || { echo "Missing value for $1" >&2; exit 2; }
      OUT_DIR="$2"; shift 2 ;;
    --done)
      [[ $# -ge 2 ]] || { echo "Missing value for $1" >&2; exit 2; }
      DONE="$2"; shift 2 ;;
    --ld-ref-dir)
      [[ $# -ge 2 ]] || { echo "Missing value for $1" >&2; exit 2; }
      LD_REF_DIR="$2"; shift 2 ;;
    --ld-ref-prefix)
      [[ $# -ge 2 ]] || { echo "Missing value for $1" >&2; exit 2; }
      LD_REF_PREFIX="$2"; shift 2 ;;
    --ld-matrix-dir)
      [[ $# -ge 2 ]] || { echo "Missing value for $1" >&2; exit 2; }
      LD_MATRIX_DIR="$2"; shift 2 ;;
    --software-dir)
      [[ $# -ge 2 ]] || { echo "Missing value for $1" >&2; exit 2; }
      SOFTWARE_DIR="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -n "$METHOD" && -n "$PRS_INPUTS_ENV" && -n "$RESOURCE_DIR" && -n "$OUT_DIR" && -n "$DONE" ]] || {
  usage >&2
  exit 2
}

[[ -f "$PRS_INPUTS_ENV" ]] || { echo "Missing PRS input env file: $PRS_INPUTS_ENV" >&2; exit 1; }

# shellcheck disable=SC1090
source "$PRS_INPUTS_ENV"

mkdir -p "$OUT_DIR" "$(dirname "$DONE")"

require_file() {
  local path="$1"
  local label="$2"
  [[ -f "$path" ]] || { echo "Missing $label: $path" >&2; exit 1; }
}

require_plink_prefix() {
  local prefix="$1"
  local label="$2"
  require_file "${prefix}.bed" "$label .bed"
  require_file "${prefix}.bim" "$label .bim"
  require_file "${prefix}.fam" "$label .fam"
}

require_sumstats_columns() {
  local path="$1"
  local label="$2"
  awk '
    NR == 1 {
      for (i = 1; i <= NF; i++) h[$i] = 1
      missing = ""
      split("SNP A1 A2 BETA SE P N", required, " ")
      for (i in required) {
        if (!(required[i] in h)) {
          missing = missing " " required[i]
        }
      }
      if (missing != "") {
        printf("Missing required columns in %s:%s\n", path, missing) > "/dev/stderr"
        exit 1
      }
      exit 0
    }
  ' path="$path" "$path" || {
    echo "Invalid $label summary statistics: $path" >&2
    exit 1
  }
}

require_file "$target_sumstats_file" "target summary statistics"
require_sumstats_columns "$target_sumstats_file" "target"
require_plink_prefix "$study_sample_plink" "ancestry-1 study sample"
[[ -z "${target_study_pheno_file:-}" ]] || require_file "$target_study_pheno_file" "ancestry-1 study phenotype"

case "$METHOD" in
  single_ct|single_prsice|single_prscs|single_ldpred2|single_lassosum2)
    ;;
  multi_ctsleb|multi_prscsx|multi_ldpred2|multi_prosper|multi_sdprs)
    require_file "$training_sumstats_file" "training summary statistics"
    require_sumstats_columns "$training_sumstats_file" "training"
    require_plink_prefix "$study_sample_plink_anc2" "ancestry-2 study sample"
    [[ -z "${training_study_pheno_file:-}" ]] || require_file "$training_study_pheno_file" "ancestry-2 study phenotype"
    ;;
  *)
    echo "Unknown PRS method: $METHOD" >&2
    exit 2
    ;;
esac

cat > "$OUT_DIR/manifest.tsv" <<EOF
key	value
method	$METHOD
resource_dir	$RESOURCE_DIR
target_sumstats_file	$target_sumstats_file
training_sumstats_file	${training_sumstats_file:-}
target_gwas_pheno_file	${target_gwas_pheno_file:-}
target_study_pheno_file	${target_study_pheno_file:-}
training_gwas_pheno_file	${training_gwas_pheno_file:-}
training_study_pheno_file	${training_study_pheno_file:-}
study_sample_plink	$study_sample_plink
study_sample_plink_anc2	${study_sample_plink_anc2:-}
reference_SNPS_bim	${reference_SNPS_bim:-}
ld_ref_dir	$LD_REF_DIR
ld_ref_prefix	$LD_REF_PREFIX
ld_matrix_dir	$LD_MATRIX_DIR
software_dir	$SOFTWARE_DIR
out_dir	$OUT_DIR
EOF

export METHOD RESOURCE_DIR OUT_DIR LD_REF_DIR LD_REF_PREFIX LD_MATRIX_DIR SOFTWARE_DIR
export target_sumstats_file training_sumstats_file study_sample_plink study_sample_plink_anc2 reference_SNPS_bim
export target_gwas_pheno_file target_study_pheno_file training_gwas_pheno_file training_study_pheno_file

if [[ -n "${PRS_METHOD_COMMAND:-}" ]]; then
  echo "Running configured command for $METHOD"
  echo "$PRS_METHOD_COMMAND" > "$OUT_DIR/command.sh"
  bash -lc "$PRS_METHOD_COMMAND"
  echo "DONE" > "$OUT_DIR/status.txt"
else
  echo "No command configured for $METHOD; inputs validated and manifest written."
  echo "READY_NO_COMMAND" > "$OUT_DIR/status.txt"
fi

touch "$DONE"
