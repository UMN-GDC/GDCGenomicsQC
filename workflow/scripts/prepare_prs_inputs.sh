#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  prepare_prs_inputs.sh --sim-dir DIR --out-dir DIR --anc1 AFR --anc2 EUR [--phenotype-index 1] [--gwas-fraction 0.5] [--gwas-n N --ld-reference-n N --tuning-n N --testing-n N] [--analysis-variant-mode none|prscsx_hm3_intersection] [--prscsx-snpinfo FILE] [--seed 42] [--plink1-bin plink] [--plink2-bin plink2] [--prs-pipeline-dir DIR] [--ld-variant-extract FILE] [--ld-prune-window-kb 500] [--ld-prune-step 50] [--ld-prune-r2 0.1] [--ld-prune-max-variants 100000]

Creates PRS pipeline inputs from simulated ancestry-specific PLINK files:
  - gwas/target_sumstats.txt and gwas/training_sumstats.txt for PRS-CSx
  - gwas/target_sumstats_singlePRS.txt and gwas/training_sumstats_singlePRS.txt for single-ancestry prs_pipeline
  - anc1_plink_files/<ANC1>_simulation_gwas.{bed,bim,fam}
  - anc1_plink_files/<ANC1>_simulation_study_sample.{bed,bim,fam}
  - anc2_plink_files/<ANC2>_simulation_gwas.{bed,bim,fam}
  - anc2_plink_files/<ANC2>_simulation_study_sample.{bed,bim,fam}
  - anc{1,2}_plink_files/<ANC>_simulation_ld_reference.{bed,bim,fam}
  - metadata/<ANC>_gwas.pheno and metadata/<ANC>_study.pheno
  - metadata/<ANC1>_tuning.{keep,pheno} and <ANC1>_testing.{keep,pheno}
    defining the shared outcome-blinded tuning/test split for every PRS method
  - anc1_plink_files/<ANC1>_simulation_study_sample_ldreduced.{bed,bim,fam}
    for LDpred2/lassosum2 (LD-pruned or matched external SNP list)
USAGE
}

SIM_DIR=""
OUT_DIR=""
ANC1="AFR"
ANC2="EUR"
PHENO_INDEX="1"
GWAS_FRACTION="0.5"
GWAS_N="0"
LD_REFERENCE_N="0"
TUNING_N="0"
TESTING_N="0"
ANALYSIS_VARIANT_MODE="none"
PRSCSX_SNPINFO=""
SEED="42"
PLINK2_BIN="plink2"
PLINK1_BIN="plink"
PRS_PIPELINE_DIR="/projects/standard/gdc/public/prs_methods/scripts/prs_pipeline"
LD_VARIANT_EXTRACT=""
LD_PRUNE_WINDOW_KB="500"
LD_PRUNE_STEP="50"
LD_PRUNE_R2="0.1"
LD_PRUNE_MAX_VARIANTS="100000"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sim-dir) SIM_DIR="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --anc1) ANC1="$2"; shift 2 ;;
    --anc2) ANC2="$2"; shift 2 ;;
    --phenotype-index) PHENO_INDEX="$2"; shift 2 ;;
    --gwas-fraction) GWAS_FRACTION="$2"; shift 2 ;;
    --gwas-n) GWAS_N="$2"; shift 2 ;;
    --ld-reference-n) LD_REFERENCE_N="$2"; shift 2 ;;
    --tuning-n) TUNING_N="$2"; shift 2 ;;
    --testing-n) TESTING_N="$2"; shift 2 ;;
    --analysis-variant-mode) ANALYSIS_VARIANT_MODE="$2"; shift 2 ;;
    --prscsx-snpinfo) PRSCSX_SNPINFO="$2"; shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;
    --plink1-bin) PLINK1_BIN="$2"; shift 2 ;;
    --plink2-bin) PLINK2_BIN="$2"; shift 2 ;;
    --prs-pipeline-dir) PRS_PIPELINE_DIR="$2"; shift 2 ;;
    --ld-variant-extract) LD_VARIANT_EXTRACT="$2"; shift 2 ;;
    --ld-prune-window-kb) LD_PRUNE_WINDOW_KB="$2"; shift 2 ;;
    --ld-prune-step) LD_PRUNE_STEP="$2"; shift 2 ;;
    --ld-prune-r2) LD_PRUNE_R2="$2"; shift 2 ;;
    --ld-prune-max-variants) LD_PRUNE_MAX_VARIANTS="$2"; shift 2 ;;
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

if ! [[ "$LD_PRUNE_MAX_VARIANTS" =~ ^[0-9]+$ ]]; then
  echo "--ld-prune-max-variants must be a non-negative integer." >&2
  exit 2
fi

for count_name in GWAS_N LD_REFERENCE_N TUNING_N TESTING_N; do
  count_value="${!count_name}"
  if ! [[ "$count_value" =~ ^[0-9]*$ ]]; then
    echo "${count_name} must be a non-negative integer." >&2
    exit 2
  fi
done

case "$ANALYSIS_VARIANT_MODE" in
  none|prscsx_hm3_intersection) ;;
  *)
    echo "--analysis-variant-mode must be none or prscsx_hm3_intersection." >&2
    exit 2
    ;;
esac

if [[ "$ANALYSIS_VARIANT_MODE" == "prscsx_hm3_intersection" ]]; then
  [[ -f "$PRSCSX_SNPINFO" ]] || {
    echo "PRS-CSx SNP metadata is required for prscsx_hm3_intersection: $PRSCSX_SNPINFO" >&2
    exit 1
  }
fi

mkdir -p "$OUT_DIR"/{gwas,anc1_plink_files,anc2_plink_files,metadata,logs,tmp}

make_split_files() {
  local fam="$1"
  local prefix="$2"
  local n
  local n_gwas n_ld_reference n_study

  n=$(awk 'NF >= 2 {n++} END {print n+0}' "$fam")
  if [[ "$n" -lt 4 ]]; then
    echo "Need at least 4 samples in $fam to create GWAS/study splits; found $n" >&2
    exit 1
  fi

  if (( GWAS_N > 0 || TUNING_N > 0 || TESTING_N > 0 )); then
    if (( GWAS_N == 0 || TUNING_N == 0 || TESTING_N == 0 )); then
      echo "--gwas-n, --tuning-n, and --testing-n must be supplied together." >&2
      exit 2
    fi
    n_gwas="$GWAS_N"
    n_ld_reference="$LD_REFERENCE_N"
    n_study=$((TUNING_N + TESTING_N))
    if (( n_gwas < 2 || n_ld_reference < 2 || TUNING_N < 2 || TESTING_N < 2 )); then
      echo "Explicit cohort sizes must each be at least two participants." >&2
      exit 2
    fi
    if (( n_gwas + n_ld_reference + n_study > n )); then
      echo "Requested GWAS/LD-reference/tuning/test cohorts exceed $prefix sample count ($n)." >&2
      exit 1
    fi
  else
    n_gwas=$(awk -v n="$n" -v f="$GWAS_FRACTION" 'BEGIN {x=int(n*f); if (x < 2) x=2; if (x > n-2) x=n-2; print x}')
    n_ld_reference=0
    n_study=$((n - n_gwas))
  fi

  awk -v seed="$SEED" 'BEGIN {srand(seed)} NF >= 2 {print rand() "\t" $1 "\t" $2}' "$fam" \
    | sort -k1,1g \
    | awk -v n_gwas="$n_gwas" -v n_ld_reference="$n_ld_reference" -v n_study="$n_study" \
        -v gwas="$OUT_DIR/metadata/${prefix}_gwas.keep" \
        -v ldref="$OUT_DIR/metadata/${prefix}_ld_reference.keep" \
        -v study="$OUT_DIR/metadata/${prefix}_study.keep" '
          NR <= n_gwas {print $2, $3 > gwas; next}
          NR <= n_gwas + n_ld_reference {print $2, $3 > ldref; next}
          NR <= n_gwas + n_ld_reference + n_study {print $2, $3 > study}'

  # Preserve legacy configurations: when no independent LD panel is requested,
  # methods receive the study cohort as their reference as before.
  if (( n_ld_reference == 0 )); then
    cp "$OUT_DIR/metadata/${prefix}_study.keep" "$OUT_DIR/metadata/${prefix}_ld_reference.keep"
  fi

  echo "${prefix} cohort sizes: GWAS=${n_gwas}, LD-reference=${n_ld_reference}, study=${n_study}"
}

make_pheno_file() {
  local subset_fam="$1"
  local simulation_fam="$2"
  local out="$3"
  local pheno_col=$((5 + PHENO_INDEX))

  # PLINK subset .fam files retain only phenotype column 1.  Recover the
  # requested replicate phenotype from the original multi-phenotype .fam.
  awk -v col="$pheno_col" '
    BEGIN {print "FID\tIID\tPHENO"}
    NR == FNR {
      if (NF >= col && $col != "NA" && $col != "-9") {
        pheno[$1 SUBSEP $2] = $col
      }
      next
    }
    {
      key = $1 SUBSEP $2
      if (key in pheno) {
        print $1 "\t" $2 "\t" pheno[key]
        n++
      } else {
        missing++
      }
    }
    END {
      if (missing > 0) {
        printf "Missing selected phenotype for %d participants in %s\\n", missing, ARGV[2] > "/dev/stderr"
        exit 1
      }
      if (n == 0) {
        printf "No phenotype values found in column %d of %s\\n", col, ARGV[1] > "/dev/stderr"
        exit 1
      }
    }
  ' "$simulation_fam" "$subset_fam" > "$out"
}

make_tuning_testing_files() {
  local study_fam="$1"
  local simulation_fam="$2"
  local anc="$3"
  local n n_tuning pheno_col
  local tuning_keep="$OUT_DIR/metadata/${anc}_tuning.keep"
  local testing_keep="$OUT_DIR/metadata/${anc}_testing.keep"
  local tuning_pheno="$OUT_DIR/metadata/${anc}_tuning.pheno"
  local testing_pheno="$OUT_DIR/metadata/${anc}_testing.pheno"

  n=$(awk 'NF >= 6 {n++} END {print n+0}' "$study_fam")
  if [[ "$n" -lt 4 ]]; then
    echo "Need at least four study participants to create tuning/testing splits; found $n" >&2
    exit 1
  fi
  if (( TUNING_N > 0 )); then
    n_tuning="$TUNING_N"
    if (( n_tuning + TESTING_N != n )); then
      echo "Study sample for $anc contains $n participants; expected tuning + testing = $((n_tuning + TESTING_N))." >&2
      exit 1
    fi
  else
    n_tuning=$((n / 2))
  fi
  pheno_col=$((5 + PHENO_INDEX))

  # The study cohort was randomly selected from the simulated population.
  # Preserve its row order because the existing CT-SLEB implementation indexes
  # its tuning and validation PRS matrices in this same order.
  awk -v n_tuning="$n_tuning" -v pheno_col="$pheno_col" \
      -v expected_n="$n" \
      -v tuning_keep="$tuning_keep" -v testing_keep="$testing_keep" \
      -v tuning_pheno="$tuning_pheno" -v testing_pheno="$testing_pheno" '
    BEGIN {
      print "FID\tIID\tPHENO" > tuning_pheno
      print "FID\tIID\tPHENO" > testing_pheno
    }
    NR == FNR {
      if (NF >= pheno_col && $pheno_col != "NA" && $pheno_col != "-9") {
        pheno[$1 SUBSEP $2] = $pheno_col
      }
      next
    }
    {
      key = $1 SUBSEP $2
      if (!(key in pheno)) {
        missing++
        next
      }
      sample_index++
    }
    sample_index <= n_tuning {
      print $1, $2 > tuning_keep
      print $1 "\t" $2 "\t" pheno[key] > tuning_pheno
      next
    }
    {
      print $1, $2 > testing_keep
      print $1 "\t" $2 "\t" pheno[key] > testing_pheno
    }
    END {
      if (missing > 0) {
        printf "Missing selected phenotype for %d study participants\\n", missing > "/dev/stderr"
        exit 1
      }
      if (sample_index != expected_n) {
        printf "Expected %d study participants with the selected phenotype; found %d\\n", expected_n, sample_index > "/dev/stderr"
        exit 1
      }
    }
  ' "$simulation_fam" "$study_fam"

  echo "Shared ${anc} tuning/testing split: ${n_tuning}/$((n - n_tuning)) participants"
}

make_analysis_variant_extracts() {
  local anc1_prefix="$1"
  local anc2_prefix="$2"
  local anc1_extract="$OUT_DIR/metadata/${ANC1}_analysis_variants.extract"
  local anc2_extract="$OUT_DIR/metadata/${ANC2}_analysis_variants.extract"
  local anc1_key_ids="$OUT_DIR/tmp/${ANC1}_prscsx_hm3_key_ids.tsv"
  local anc2_key_ids="$OUT_DIR/tmp/${ANC2}_prscsx_hm3_key_ids.tsv"

  if [[ "$ANALYSIS_VARIANT_MODE" == "none" ]]; then
    awk 'NF {print $2}' "${anc1_prefix}.bim" > "$anc1_extract"
    awk 'NF {print $2}' "${anc2_prefix}.bim" > "$anc2_extract"
    return 0
  fi

  # Both ancestry panels retain their native IDs at loci represented in the
  # PRS-CSx 1KG HapMap3 reference, yielding a common analysis variant set.
  awk 'NR == FNR {if (FNR > 1) ref[$1 ":" $3] = 1; next} {key = $1 ":" $4; if (key in ref) print key "\t" $2}' \
    "$PRSCSX_SNPINFO" "${anc1_prefix}.bim" > "$anc1_key_ids"
  awk 'NR == FNR {if (FNR > 1) ref[$1 ":" $3] = 1; next} {key = $1 ":" $4; if (key in ref) print key "\t" $2}' \
    "$PRSCSX_SNPINFO" "${anc2_prefix}.bim" > "$anc2_key_ids"
  awk -v anc1_out="$anc1_extract" -v anc2_out="$anc2_extract" '
    NR == FNR {anc1[$1] = $2; next}
    $1 in anc1 {print anc1[$1] > anc1_out; print $2 > anc2_out}
  ' "$anc1_key_ids" "$anc2_key_ids"

  local n_anc1 n_anc2
  n_anc1=$(awk 'NF {n++} END {print n+0}' "$anc1_extract")
  n_anc2=$(awk 'NF {n++} END {print n+0}' "$anc2_extract")
  if (( n_anc1 < 100000 || n_anc1 != n_anc2 )); then
    echo "Shared PRS-CSx/HM3 variant extraction failed: ${ANC1}=${n_anc1}, ${ANC2}=${n_anc2}" >&2
    exit 1
  fi
  echo "Shared PRS-CSx/HM3 positional intersection: ${n_anc1} variants per ancestry"
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
  local ld_reference_prefix="$OUT_DIR/${subdir}/${anc}_simulation_ld_reference"
  local gwas_pheno="$OUT_DIR/metadata/${anc}_gwas.pheno"
  local study_pheno="$OUT_DIR/metadata/${anc}_study.pheno"
  local glm_prefix="$OUT_DIR/gwas/${anc}_glm"

  for ext in bed bim fam; do
    [[ -f "${sim_prefix}.${ext}" ]] || { echo "Missing ${sim_prefix}.${ext}" >&2; exit 1; }
  done

  "$PLINK2_BIN_RESOLVED" --bfile "$sim_prefix" \
    --extract "$OUT_DIR/metadata/${anc}_analysis_variants.extract" \
    --make-bed --out "$out_prefix"

  make_split_files "${out_prefix}.fam" "$anc"

  "$PLINK2_BIN_RESOLVED" --bfile "$out_prefix" --keep "$OUT_DIR/metadata/${anc}_gwas.keep" --make-bed --out "$gwas_prefix"
  "$PLINK2_BIN_RESOLVED" --bfile "$out_prefix" --keep "$OUT_DIR/metadata/${anc}_ld_reference.keep" --make-bed --out "$ld_reference_prefix"
  "$PLINK2_BIN_RESOLVED" --bfile "$out_prefix" --keep "$OUT_DIR/metadata/${anc}_study.keep" --make-bed --out "$study_prefix"

  make_pheno_file "${gwas_prefix}.fam" "${sim_prefix}.fam" "$gwas_pheno"
  make_pheno_file "${study_prefix}.fam" "${sim_prefix}.fam" "$study_pheno"
  run_glm "$gwas_prefix" "$gwas_pheno" "$glm_prefix"
}

make_ld_reduced_study_set() {
  local study_prefix="$1"
  local anc="$2"
  local reduced_prefix="$OUT_DIR/anc1_plink_files/${anc}_simulation_study_sample_ldreduced"
  local prune_prefix="$OUT_DIR/tmp/${anc}_ldprune"
  local selected_ids="$OUT_DIR/metadata/${anc}_ld_variant_ids.txt"
  local n_selected

  if [[ -n "$LD_VARIANT_EXTRACT" ]]; then
    [[ -f "$LD_VARIANT_EXTRACT" ]] || {
      echo "LD variant extract file not found: $LD_VARIANT_EXTRACT" >&2
      exit 1
    }
    awk 'NF {print $1}' "$LD_VARIANT_EXTRACT" > "$selected_ids"
  else
    "$PLINK2_BIN_RESOLVED" \
      --bfile "$study_prefix" \
      --indep-pairwise "$LD_PRUNE_WINDOW_KB" "$LD_PRUNE_STEP" "$LD_PRUNE_R2" \
      --out "$prune_prefix"

    if (( LD_PRUNE_MAX_VARIANTS > 0 )); then
      # Select a reproducible genome-wide subset instead of favoring early chromosomes.
      awk -v seed="$SEED" 'BEGIN {srand(seed)} {print rand() "\t" $1}' "${prune_prefix}.prune.in" \
        | sort -k1,1g \
        | awk -v limit="$LD_PRUNE_MAX_VARIANTS" 'NR <= limit {print $2}' > "$selected_ids"
    else
      cp "${prune_prefix}.prune.in" "$selected_ids"
    fi
  fi

  n_selected=$(awk 'NF {n++} END {print n+0}' "$selected_ids")
  if [[ "$n_selected" -lt 2 ]]; then
    echo "LD-reduced set contains fewer than two variants: $selected_ids" >&2
    exit 1
  fi

  "$PLINK2_BIN_RESOLVED" \
    --bfile "$study_prefix" \
    --extract "$selected_ids" \
    --make-bed \
    --out "$reduced_prefix"

  echo "LD-reduced study set contains $n_selected variants: $reduced_prefix"
}

make_analysis_variant_extracts "$SIM_DIR/${ANC1}_simulation" "$SIM_DIR/${ANC2}_simulation"
process_ancestry "$ANC1" "anc1_plink_files"
make_tuning_testing_files \
  "$OUT_DIR/anc1_plink_files/${ANC1}_simulation_study_sample.fam" \
  "$SIM_DIR/${ANC1}_simulation.fam" \
  "$ANC1"
make_ld_reduced_study_set "$OUT_DIR/anc1_plink_files/${ANC1}_simulation_study_sample" "$ANC1"
format_sumstats "$OUT_DIR/gwas/${ANC1}_glm" "$OUT_DIR/gwas/target_sumstats.txt"
format_single_prs_sumstats "$OUT_DIR/gwas/target_sumstats.txt" "$OUT_DIR/gwas/target_sumstats_singlePRS.txt"

process_ancestry "$ANC2" "anc2_plink_files"
format_sumstats "$OUT_DIR/gwas/${ANC2}_glm" "$OUT_DIR/gwas/training_sumstats.txt"
format_single_prs_sumstats "$OUT_DIR/gwas/training_sumstats.txt" "$OUT_DIR/gwas/training_sumstats_singlePRS.txt"

cohort_manifest="$OUT_DIR/metadata/cohort_allocation.tsv"
{
  printf 'ancestry\tcohort\tn\n'
  for anc in "$ANC1" "$ANC2"; do
    for cohort in gwas ld_reference study; do
      printf '%s\t%s\t%s\n' "$anc" "$cohort" "$(awk 'END {print NR+0}' "$OUT_DIR/metadata/${anc}_${cohort}.keep")"
    done
  done
  for cohort in tuning testing; do
    printf '%s\t%s\t%s\n' "$ANC1" "$cohort" "$(awk 'END {print NR+0}' "$OUT_DIR/metadata/${ANC1}_${cohort}.keep")"
  done
  printf '%s\tshared_prscsx_hm3_variants\t%s\n' "$ANC1" "$(awk 'END {print NR+0}' "$OUT_DIR/metadata/${ANC1}_analysis_variants.extract")"
} > "$cohort_manifest"

cat > "$OUT_DIR/prs_inputs.env" <<EOF
path_data_root="$OUT_DIR"
target_sumstats_file="$OUT_DIR/gwas/target_sumstats.txt"
training_sumstats_file="$OUT_DIR/gwas/training_sumstats.txt"
target_single_prs_sumstats_file="$OUT_DIR/gwas/target_sumstats_singlePRS.txt"
training_single_prs_sumstats_file="$OUT_DIR/gwas/training_sumstats_singlePRS.txt"
target_gwas_pheno_file="$OUT_DIR/metadata/${ANC1}_gwas.pheno"
target_study_pheno_file="$OUT_DIR/metadata/${ANC1}_study.pheno"
target_tuning_keep_file="$OUT_DIR/metadata/${ANC1}_tuning.keep"
target_testing_keep_file="$OUT_DIR/metadata/${ANC1}_testing.keep"
target_tuning_pheno_file="$OUT_DIR/metadata/${ANC1}_tuning.pheno"
target_testing_pheno_file="$OUT_DIR/metadata/${ANC1}_testing.pheno"
training_gwas_pheno_file="$OUT_DIR/metadata/${ANC2}_gwas.pheno"
training_study_pheno_file="$OUT_DIR/metadata/${ANC2}_study.pheno"
reference_SNPS_bim="$OUT_DIR/anc1_plink_files/${ANC1}_simulation_study_sample"
study_sample_plink="$OUT_DIR/anc1_plink_files/${ANC1}_simulation_study_sample"
target_ld_reference_plink="$OUT_DIR/anc1_plink_files/${ANC1}_simulation_ld_reference"
aux_ld_reference_plink="$OUT_DIR/anc2_plink_files/${ANC2}_simulation_ld_reference"
cohort_allocation_file="$cohort_manifest"
ld_study_sample="$OUT_DIR/anc1_plink_files/${ANC1}_simulation_study_sample_ldreduced"
study_sample_plink_anc2="$OUT_DIR/anc2_plink_files/${ANC2}_simulation_study_sample"
EOF

N_GWAS_ANC1=$(awk 'END {print NR}' "$OUT_DIR/metadata/${ANC1}_gwas.keep")

cat > "$OUT_DIR/prs_prscsx_generated.conf" <<EOF
path_code="/projects/standard/gdc/public/prs_methods/scripts/PRScsx"
path_data_root="$OUT_DIR"
path_ref_dir="/projects/standard/gdc/public/prs_methods/ref/ref_PRScsx/1kg_ref"
path_plink2="$PLINK2_BIN_RESOLVED"

anc1="$ANC1"
anc2="$ANC2"

target_sumstats_file="\${path_data_root}/gwas/target_sumstats.txt"
training_sumstats_file="\${path_data_root}/gwas/training_sumstats.txt"
target_study_pheno_file="\${path_data_root}/metadata/${ANC1}_study.pheno"
target_tuning_keep_file="\${path_data_root}/metadata/${ANC1}_tuning.keep"
target_testing_keep_file="\${path_data_root}/metadata/${ANC1}_testing.keep"
target_tuning_pheno_file="\${path_data_root}/metadata/${ANC1}_tuning.pheno"
target_testing_pheno_file="\${path_data_root}/metadata/${ANC1}_testing.pheno"
training_study_pheno_file="\${path_data_root}/metadata/${ANC2}_study.pheno"
output_dir="$OUT_DIR"

reference_SNPS_bim="\${path_data_root}/anc1_plink_files/${ANC1}_simulation_study_sample"
study_sample_plink="\${path_data_root}/anc1_plink_files/${ANC1}_simulation_study_sample"
target_ld_reference_plink="\${path_data_root}/anc1_plink_files/${ANC1}_simulation_ld_reference"
aux_ld_reference_plink="\${path_data_root}/anc2_plink_files/${ANC2}_simulation_ld_reference"
study_sample_plink_anc2="\${path_data_root}/anc2_plink_files/${ANC2}_simulation_study_sample"

prs_pipeline="$PRS_PIPELINE_DIR"
EOF

cat > "$OUT_DIR/prs_single_ancestry_${ANC1}_generated.conf" <<EOF
path_data="$OUT_DIR"
path_repo="$PRS_PIPELINE_DIR"
path_plink1="$PLINK1_BIN"
path_plink2="$PLINK2_BIN_RESOLVED"

summary_stats_file="\${path_data}/gwas/target_sumstats_singlePRS.txt"
bim_file_path="\${path_data}/anc1_plink_files/${ANC1}_simulation_study_sample.bim"
study_sample="\${path_data}/anc1_plink_files/${ANC1}_simulation_study_sample"
clump_reference_sample="\${path_data}/anc1_plink_files/${ANC1}_simulation_ld_reference"
ld_study_sample="\${path_data}/anc1_plink_files/${ANC1}_simulation_study_sample_ldreduced"
output_path="\${path_data}/single_ancestry_${ANC1}"

# C+T derives PCs from the matching study sample when this is unset.
gwas_pca_eigenvec_file=""
tuning_keep_file="\${path_data}/metadata/${ANC1}_tuning.keep"
testing_keep_file="\${path_data}/metadata/${ANC1}_testing.keep"
tuning_pheno_file="\${path_data}/metadata/${ANC1}_tuning.pheno"
testing_pheno_file="\${path_data}/metadata/${ANC1}_testing.pheno"

n_total_gwas=$N_GWAS_ANC1
cohort_allocation_file="\${path_data}/metadata/cohort_allocation.tsv"
EOF

echo "PRS inputs written under $OUT_DIR"
echo "PRS-CSx config: $OUT_DIR/prs_prscsx_generated.conf"
echo "Single-ancestry config: $OUT_DIR/prs_single_ancestry_${ANC1}_generated.conf"
