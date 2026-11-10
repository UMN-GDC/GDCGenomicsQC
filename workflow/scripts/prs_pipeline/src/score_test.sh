#!/bin/bash
#
# Score test data using trained PRS models and evaluate performance
# Usage: score_test.sh [options]

set -eu

usage() {
  cat <<EOF
Usage: $0 [options]

Required:
  --test-bfile <prefix>     PLINK prefix for test data
  --train-out-dir <dir>     Pipeline methods output directory (e.g., output_path/prs_pipeline)
  --sumstats <file>         GWAS summary stats file (needed for C+T)
  --path-repo <dir>         Path to prs_pipeline repo

Optional:
  --pheno-file <file>       Phenotype file (FID IID phenotype). Default: .fam column 6
  --test-pca-file <file>    PCA eigenvec file for test samples. Default: PLINK --pca 6
  --binary-flag <T/F>       Binary phenotype flag for PRSice2 (default: F)
  --ran-ct <true/false>     Whether C+T was run (default: false)
  --ran-ldpred2 <true/false> Whether LDpred2 was run (default: false)
  --ran-lassosum2 <true/false> Whether lassosum2 was run (default: false)
  --ran-prsice2 <true/false> Whether PRSice2 was run (default: false)
  --ran-prscsx <true/false> Whether PRS-CSx was run (default: false)
  --ran-viprs <true/false> Whether VIPRS was run (default: false)
  --ran-prosper <true/false> Whether PROSPER was run (default: false)
  -h, --help                Show this help
EOF
  exit 1
}

# --- Parse args ---
test_bfile=""
train_out_dir=""
pheno_file=""
test_pca_file=""
sumstats=""
path_repo=""
binary_flag=F
ran_ct=false
ran_ldpred2=false
ran_lassosum2=false
ran_prsice2=false
ran_prscsx=false
ran_viprs=false
ran_prosper=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --test-bfile)     test_bfile="$2";     shift 2 ;;
    --train-out-dir)  train_out_dir="$2";  shift 2 ;;
    --pheno-file)     pheno_file="$2";     shift 2 ;;
    --test-pca-file)  test_pca_file="$2";  shift 2 ;;
    --sumstats)       sumstats="$2";       shift 2 ;;
    --path-repo)      path_repo="$2";      shift 2 ;;
    --binary-flag)    binary_flag="$2";    shift 2 ;;
    --ran-ct)         ran_ct="$2";         shift 2 ;;
    --ran-ldpred2)    ran_ldpred2="$2";    shift 2 ;;
    --ran-lassosum2)  ran_lassosum2="$2";  shift 2 ;;
    --ran-prsice2)    ran_prsice2="$2";    shift 2 ;;
    --ran-prscsx)     ran_prscsx="$2";     shift 2 ;;
    --ran-viprs)      ran_viprs="$2";      shift 2 ;;
    --ran-prosper)    ran_prosper="$2";    shift 2 ;;
    -h|--help)        usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

if [[ -z "$test_bfile" || -z "$train_out_dir" || -z "$sumstats" || -z "$path_repo" ]]; then
  echo "ERROR: --test-bfile, --train-out-dir, --sumstats, and --path-repo are required."
  usage
fi

eval_dir="${train_out_dir}/test_evaluation"
mkdir -p "$eval_dir"

echo "[score_test] Evaluating PRS models on test data: $test_bfile"
echo "[score_test] Output dir: $eval_dir"

# --- 1. PCA for test samples ---
if [[ -n "$test_pca_file" ]]; then
  echo "[score_test] Using user-provided PCA file: $test_pca_file"
  pca_file="$test_pca_file"
else
  echo "[score_test] Computing PCA on test sample..."
  plink --bfile "$test_bfile" --pca 6 --allow-no-sex --out "${eval_dir}/pca"
  pca_file="${eval_dir}/pca.eigenvec"
fi

# --- 2. Phenotype for test samples ---
align_pheno() {
  local src="$1"
  local ref_fam="${test_bfile}.fam"
  local out="${eval_dir}/pheno.txt"
  # Detect column count: 3 = "FID IID pheno", 2 = "ID pheno"
  local ncols
  ncols=$(awk 'NR==1 {print NF}' "$src")
  echo "[score_test]   Detected ${ncols}-column phenotype file"
  # Build lookup by IID; handles both 2-col (ID pheno) and 3+col (FID IID pheno...) formats
  if [[ "$ncols" -ge 3 ]]; then
    awk 'NR==FNR {
         if (FNR==1 && ($1 ~ /^[Ff][Ii][Dd]/ || $2 ~ /[A-Za-z]/)) next
         pheno[$2]=$3; next
       }
       FNR==1 {print "FID\tIID\tphenotype"; next}
       $2 in pheno {print $1, $2, pheno[$2]}' \
      OFS="\t" "$src" "$ref_fam" > "$out"
  else
    awk 'NR==FNR {
         if (FNR==1 && ($1 ~ /^[Ff][Ii][Dd]/ || $1 ~ /[A-Za-z]/)) next
         pheno[$1]=$2; next
       }
       FNR==1 {print "FID\tIID\tphenotype"; next}
       $2 in pheno {print $1, $2, pheno[$2]}' \
      OFS="\t" "$src" "$ref_fam" > "$out"
  fi
  local n
  n=$(( $(wc -l < "$out") - 1 ))
  echo "[score_test]   Aligned ${n} samples"
  # Fallback: if no IIDs matched, use .fam column 6
  if [[ "$n" -eq 0 ]]; then
    echo "[score_test]   WARNING: No IID matches — falling back to .fam column 6"
    awk 'BEGIN{print "FID\tIID\tphenotype"} {print $1, $2, $6}' OFS="\t" "$ref_fam" > "$out"
    echo "[score_test]   Using $(( $(wc -l < "$out") - 1 )) samples from .fam"
  fi
}

if [[ -n "$pheno_file" ]]; then
  echo "[score_test] Using user-provided phenotype file: $pheno_file"
  align_pheno "$pheno_file"
  pheno_input="${eval_dir}/pheno.txt"
else
  echo "[score_test] Extracting phenotype from .fam column 6"
  awk 'BEGIN{print "FID\tIID\tphenotype"} {print $1, $2, $6}' OFS="\t" "${test_bfile}.fam" > "${eval_dir}/pheno.txt"
  pheno_input="${eval_dir}/pheno.txt"
fi

# --- 3. Score + evaluate each trained method ---
evaluate() {
  local label="$1"
  local out_prefix="$2"
  local profile_file="${eval_dir}/${label}.profile"

  if [[ ! -f "$profile_file" ]]; then
    echo "[score_test] WARNING: Profile not found: $profile_file — skipping ${label}"
    return
  fi

  Rscript "${path_repo}/src/evaluate_test.R" \
    "$profile_file" \
    "$pheno_input" \
    "$pca_file" \
    "${eval_dir}/${label}_results.txt" \
    "${eval_dir}/${label}_scores.txt"
}

# --- C+T ---
if [[ "$ran_ct" == true ]]; then
  ct_dir="${train_out_dir}/CT"
  ct_results="${ct_dir}/CT_prs_results.txt"
  if [[ -f "$ct_results" ]]; then
    echo "[score_test] Evaluating C+T..."
    # Find best threshold (max R2)
    best_p=$(awk 'NR>1 && $2+0 > max {max=$2; best=$1} END {print best}' "$ct_results")
    if [[ -z "$best_p" ]]; then
      echo "[score_test] WARNING: Could not parse best p-value from $ct_results — skipping C+T"
    else
      echo "[score_test]   Best C+T p-value threshold: $best_p"
      # Create single-range file at best threshold
      echo "${best_p} 0 ${best_p}" > "${eval_dir}/CT_best_range.txt"
      # Recreate SNP.pvalue from sumstats
      awk 'NR==1 {print "SNP pvalue"; next} {print $1, $8}' OFS="\t" "$sumstats" > "${eval_dir}/CT_SNP.pvalue"
      # Extract clumped SNP list from training temp dir
      valid_snp="${ct_dir}/temp/temp.valid.snp"
      if [[ ! -f "$valid_snp" ]]; then
        echo "[score_test] WARNING: Clumped SNP list not found: $valid_snp — skipping C+T"
      else
        plink --bfile "$test_bfile" \
          --score "$sumstats" 1 4 6 header \
          --q-score-range "${eval_dir}/CT_best_range.txt" "${eval_dir}/CT_SNP.pvalue" \
          --extract "$valid_snp" \
          --allow-no-sex \
          --out "${eval_dir}/CT"
        # q-score-range names the profile CT.{range_name}.profile; normalize to CT.profile
        actual_ct_profile=$(ls "${eval_dir}"/CT.*.profile 2>/dev/null | head -1)
        if [[ -f "$actual_ct_profile" ]]; then
          cp "$actual_ct_profile" "${eval_dir}/CT.profile"
        fi
        evaluate "CT" "${eval_dir}/CT"
      fi
    fi
  else
    echo "[score_test] WARNING: C+T results not found: $ct_results — skipping"
  fi
fi

# --- LDpred2 ---
if [[ "$ran_ldpred2" == true ]]; then
  ldpred2_base="${train_out_dir}/LDpred2/prs_method"
  echo "[score_test] Evaluating LDpred2..."
  for model in inf grid; do
    weights="${ldpred2_base}_${model}_weights.txt"
    if [[ -f "$weights" ]]; then
      plink --bfile "$test_bfile" \
        --score "$weights" 1 2 3 header \
        --allow-no-sex \
        --out "${eval_dir}/LDpred2_${model}"
      evaluate "LDpred2_${model}" "${eval_dir}/LDpred2_${model}"
    else
      echo "[score_test] WARNING: LDpred2 ${model} weights not found: $weights — skipping"
    fi
  done
fi

# --- lassosum2 ---
if [[ "$ran_lassosum2" == true ]]; then
  lassosum2_base="${train_out_dir}/lassosum2/prs_method"
  echo "[score_test] Evaluating lassosum2..."
  weights="${lassosum2_base}_weights.txt"
  if [[ -f "$weights" ]]; then
    plink --bfile "$test_bfile" \
      --score "$weights" 1 2 3 header \
      --allow-no-sex \
      --out "${eval_dir}/lassosum2"
    evaluate "lassosum2" "${eval_dir}/lassosum2"
  else
    echo "[score_test] WARNING: lassosum2 weights not found: $weights — skipping"
  fi
fi

# --- PRSice2 ---
if [[ "$ran_prsice2" == true ]]; then
  echo "[score_test] Evaluating PRSice2 on test data..."
  prsice2_base="${train_out_dir}/PRSice2/prs_method"
  prsice2_prsice="${prsice2_base}/PRSice2_outputs.prsice"
  prsice2_snps="${prsice2_base}/PRSice2_outputs.snps"
  if [[ ! -f "$prsice2_prsice" ]]; then
    echo "[score_test] WARNING: PRSice2 results not found: $prsice2_prsice — skipping PRSice2"
  elif [[ ! -f "$prsice2_snps" ]]; then
    echo "[score_test] WARNING: PRSice2 SNP list not found: $prsice2_snps — skipping PRSice2"
  else
    # Diagnose .prsice content
    prsice2_nlines=$(wc -l < "$prsice2_prsice")
    prsice2_header=$(head -1 "$prsice2_prsice")
    prsice2_ncols=$(awk '{print NF; exit}' "$prsice2_prsice")
    echo "[score_test]   .prsice: ${prsice2_nlines} lines, ${prsice2_ncols} cols, header: ${prsice2_header}"
    # Try .prsice (max R2) and .summary (single best row) as fallback
    best_p=$(Rscript --vanilla -e '
      f <- commandArgs(trailingOnly = TRUE)[1]
      d <- tryCatch(read.table(f, header = TRUE), error = function(e) NULL)
      if (!is.null(d) && nrow(d) > 0) {
        r2 <- grep("^R2$", colnames(d), ignore.case = TRUE)[1]
        th <- grep("Threshold", colnames(d), ignore.case = TRUE)[1]
        if (!is.na(r2) && !is.na(th))
          cat(as.character(d[which.max(d[[r2]]), th]))
      }
    ' "$prsice2_prsice")
    if [[ -z "$best_p" || ! "$best_p" =~ ^[0-9] ]]; then
      prsice2_summary="${prsice2_base}/PRSice2_outputs.summary"
      echo "[score_test]   .prsice gave no threshold — trying .summary"
      if [[ -f "$prsice2_summary" ]]; then
        best_p=$(Rscript --vanilla -e '
          f <- commandArgs(trailingOnly = TRUE)[1]
          d <- tryCatch(read.table(f, header = TRUE), error = function(e) NULL)
          if (!is.null(d) && nrow(d) > 0) {
            th <- grep("Threshold", colnames(d), ignore.case = TRUE)[1]
            if (!is.na(th)) cat(as.character(d[1, th]))
          }
        ' "$prsice2_summary")
      fi
    fi
    if [[ -z "$best_p" || ! "$best_p" =~ ^[0-9] ]]; then
      echo "[score_test] WARNING: Could not parse best p-value from .prsice or .summary — skipping PRSice2"
    else
      echo "[score_test]   Best PRSice2 p-value threshold from training: $best_p"
      # Score test data using PRSice2 with training-derived parameters (per official docs)
      # --no-clump: use pre-clumped SNPs from training
      # --extract: restrict to training SNP list
      # --fastscore + --bar-levels: PRS only at training's best threshold
      # --no-regress: output PRS without regression on test data
      prsice2_test_out="${eval_dir}/PRSice2"
      mkdir -p "$prsice2_test_out"
      PRSice \
        --base "$sumstats" \
        --target "$test_bfile" \
        --binary-target "${binary_flag:-F}" \
        --beta \
        --stat beta \
        --score avg \
        --no-clump \
        --extract "$prsice2_snps" \
        --bar-levels "$best_p" \
        --fastscore \
        --no-regress \
        --out "${prsice2_test_out}/PRSice2_test" 2>/dev/null
      if [[ -f "${prsice2_test_out}/PRSice2_test.best" ]]; then
        awk 'NR==1 {print "FID\tIID\tSCORE"} NR>1 {print $1"\t"$2"\t"$3}' \
          "${prsice2_test_out}/PRSice2_test.best" > "${eval_dir}/PRSice2.profile"
        evaluate "PRSice2" "${eval_dir}/PRSice2"
      else
        echo "[score_test] WARNING: PRSice2 .best not found — falling back to PLINK scoring"
        echo "${best_p} 0 ${best_p}" > "${eval_dir}/PRSice2_best_range.txt"
        awk 'NR==1 {print "SNP pvalue"; next} {print $1, $8}' OFS="\t" "$sumstats" > "${eval_dir}/PRSice2_SNP.pvalue"
        plink --bfile "$test_bfile" \
          --score "$sumstats" 1 4 6 header \
          --q-score-range "${eval_dir}/PRSice2_best_range.txt" "${eval_dir}/PRSice2_SNP.pvalue" \
          --extract "$prsice2_snps" \
          --allow-no-sex \
          --out "${eval_dir}/PRSice2"
        actual_prsice2_profile=$(ls "${eval_dir}"/PRSice2.*.profile 2>/dev/null | head -1)
        if [[ -f "$actual_prsice2_profile" ]]; then
          cp "$actual_prsice2_profile" "${eval_dir}/PRSice2.profile"
        fi
        evaluate "PRSice2" "${eval_dir}/PRSice2"
      fi
    fi
  fi
fi

# --- PRS-CSx ---
if [[ "$ran_prscsx" == true ]]; then
  echo "[score_test] Evaluating PRS-CSx on test data..."
  prscsx_weights=("${train_out_dir}"/PRScsx/*_combined_weights.txt)
  if [[ -e "${prscsx_weights[0]}" ]]; then
    for weights in "${prscsx_weights[@]}"; do
      anc_label=$(basename "$weights" "_combined_weights.txt")
      label="PRScsx_${anc_label}"
      echo "[score_test]   Scoring test data with ${anc_label} weights"
      # Strip the (3-col) header line — data rows have 6 cols (CHR SNP POS A1 A2 BETA)
      awk 'NR>1' "$weights" > "${eval_dir}/${label}.weights"
      plink --bfile "$test_bfile" \
        --score "${eval_dir}/${label}.weights" 2 4 6 \
        --allow-no-sex \
        --out "${eval_dir}/${label}"
      evaluate "$label" "${eval_dir}/${label}"
    done
  else
    echo "[score_test] WARNING: No PRS-CSx combined weights found in ${train_out_dir}/PRScsx — skipping"
  fi
fi

# --- VIPRS ---
if [[ "$ran_viprs" == true ]]; then
  echo "[score_test] Evaluating VIPRS on test data..."
  viprs_fit_file="${train_out_dir}/viprs/VIPRS_EM.fit.gz"
  if [[ ! -f "$viprs_fit_file" ]]; then
    echo "[score_test] WARNING: VIPRS fit file not found: $viprs_fit_file — skipping"
  else
    viprs_score -f "$viprs_fit_file" \
      --bfile "$test_bfile" \
      --output-file "${eval_dir}/VIPRS_test" 2>/dev/null
    if [[ ! -f "${eval_dir}/VIPRS_test.prs" ]]; then
      echo "[score_test] WARNING: VIPRS scoring produced no .prs file — skipping"
    else
      # Convert VIPRS .prs (FID IID PRS) to profile format (FID IID SCORE),
      # recovering FID from the test .fam
      echo -e "FID\tIID\tSCORE" > "${eval_dir}/VIPRS.profile"
      awk 'NR==FNR { fid[$2]=$1; next }
           FNR==1 && $1 !~ /^[0-9]/ { next }
           { print ($2 in fid ? fid[$2] : $1), $2, $3 }' \
        OFS="\t" "${test_bfile}.fam" "${eval_dir}/VIPRS_test.prs" >> "${eval_dir}/VIPRS.profile"
      evaluate "VIPRS" "${eval_dir}/VIPRS"
    fi
  fi
fi

# --- PROSPER ---
if [[ "$ran_prosper" == true ]]; then
  echo "[score_test] Evaluating PROSPER on test data..."
  prosper_prs=$(find "${train_out_dir}/PROSPER" -name "PROSPER_prs_file.txt" -path "*/after_ensemble_*" 2>/dev/null | head -1)
  if [[ -z "$prosper_prs" || ! -f "$prosper_prs" ]]; then
    echo "[score_test] WARNING: PROSPER PRS weights not found in ${train_out_dir}/PROSPER/after_ensemble_*/PROSPER_prs_file.txt — skipping"
  else
    # PROSPER_prs_file.txt columns: rsid, a1, a0, weight (with header)
    # Convert to plink --score format: SNP A1 weight
    awk 'NR==1 {print "SNP", "A1", "weight"; next} {print $1, $2, $4}' OFS="\t" "$prosper_prs" \
      > "${eval_dir}/PROSPER.weights"
    n_weights=$(wc -l < "${eval_dir}/PROSPER.weights")
    if [[ "$n_weights" -le 1 ]]; then
      echo "[score_test] WARNING: PROSPER PRS file has no SNP weights (${prosper_prs}) — skipping"
    else
      echo "[score_test]   Using PROSPER weights: $prosper_prs (${n_weights} SNPs)"
      plink --bfile "$test_bfile" \
        --score "${eval_dir}/PROSPER.weights" 1 2 3 header \
        --allow-no-sex \
        --out "${eval_dir}/PROSPER"
      evaluate "PROSPER" "${eval_dir}/PROSPER"
    fi
  fi
fi

echo "[score_test] Done. Results in: $eval_dir"
