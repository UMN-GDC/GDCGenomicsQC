#!/usr/bin/env bash
#
# smoke_single_ancestry.sh
#
# End-to-end SLURM smoke test for the five single-ancestry PRS methods
# (CT, PRSice2, PRS-CS, LDpred2, lassosum2) plus the joint multi-ancestry
# PRS-CSx method (multi_prscsx) in GDCGenomicsQC.
#
# What it does:
#   1. Preflights tools (snakemake, Rscript, plink/plink2, PRS-CS python).
#   2. Fabricates two small single-ancestry toy PLINK1 sets (AFR + EUR) whose
#      variants are REAL chr22 HapMap3 SNPs (so PRS-CS LD-overlap works).
#   3. Seeds the phenotypeSimulation outputs (bypasses the QC/sim QC graph).
#   4. Writes a full config, dry-runs the graph, then RUNS
#      `run_singleAncestryPRSPipelines` and `runMultiAncestryPRSCSx` locally
#      within this allocation.
#   5. Verifies all five single method outputs + `.done` markers, the
#      multi_prscsx outputs for both ancestries, and (Phase 3) prints a table.
#
# Run via:
#   sbatch workflow/scripts/smoke_single_ancestry.sh
# or directly (no SLURM):
#   bash workflow/scripts/smoke_single_ancestry.sh
#
# Overridable env vars (all optional):
#   GDCQC_SMOKE_WORK           workspace root (default: /scratch.global/baron063/testing/GDCQC_PRS_integration)
#   GDCQC_SMOKE_KEEP=1         reuse an existing workspace instead of wiping it
#   GDCQC_SMOKE_DRYRUN=1       config-validating dry-run only, no execution
#   GDCQC_SMOKE_SNPS           variants per toy set (default 300)
#   GDCQC_SMOKE_SAMPLES        samples per toy set (default 100)
#   GDCQC_SMOKE_SEED           RNG seed (default 42)
#   GDCQC_SMOKE_TRACE=1        add --printshellcmds to the snakemake runs
#
#SBATCH --job-name=gdcqc-prs-smoke
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=gdcqc_prs_smoke_%j.log
#SBATCH --error=gdcqc_prs_smoke_%j.log

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
WORKFLOW="$REPO_ROOT/workflow"
GENERATOR="$SCRIPT_DIR/make_smoke_toy.py"

START_EPOCH="$(date +%s)"
log()  { printf '\n[%s] %s\n' "$(date '+%F %T')" "$*"; }
section() {
    echo
    echo "============================================================"
    echo "  $*"
    echo "============================================================"
}
die() { printf 'ERROR: %s\n' "$1" >&2; exit 1; }

runas() {  # runas <label> <cmd...>
    local label="$1"; shift
    log "RUN [$label]: $*"
    "$@"
    local rc=$?
    [[ $rc -eq 0 ]] || die "[$label] failed with rc=$rc (see output above)"
}

# ---------------------------------------------------------------- workspace
if [[ ${GDCQC_SMOKE_WORK:-} != "" ]]; then
    WORK_BASE="$GDCQC_SMOKE_WORK"
elif [[ $USER == baron063 || -w /scratch.global/baron063 ]]; then
    WORK_BASE="/scratch.global/baron063/testing/GDCQC_PRS_integration"
else
    WORK_BASE="$HOME/gdcqc_smoke"
fi
OUT_DIR="$WORK_BASE/OUT"
SIM_INPUTS="$WORK_BASE/sim_inputs"
RESOURCES_DIR="$WORK_BASE/prs_resources"
CONFIG="$WORK_BASE/smoke_config.yaml"
RUN_LOG="$WORK_BASE/snakemake_run.log"
ANC1=AFR
ANC2=EUR
PRS_OUT="$OUT_DIR/prs_inputs/${ANC1}_${ANC2}"
METHOD_RUN="$PRS_OUT/method_runs"
NSNPS="${GDCQC_SMOKE_SNPS:-300}"
NSAMPLES="${GDCQC_SMOKE_SAMPLES:-100}"
SEED="${GDCQC_SMOKE_SEED:-42}"
CPUS="${GDCQC_SMOKE_CPUS:-${SLURM_CPUS_ON_NODE:-8}}"
TRACE_ARGS=""
[[ ${GDCQC_SMOKE_TRACE:-0} == 1 ]] && TRACE_ARGS="--printshellcmds"

section "GDCGenomicsQC single-ancestry PRS smoke"
log "repo root        : $REPO_ROOT"
log "workflow dir     : $WORKFLOW"
log "workspace        : $WORK_BASE"
log "toy SNPs/samples : $NSNPS / $NSAMPLES x2 (seed $SEED)"
log "cores            : $CPUS"
log "started          : $(date)"

for f in "$WORKFLOW/Snakefile" "$GENERATOR" "$WORKFLOW/scripts/prs_pipeline/run_single_ancestry_PRS_pipeline.sh"; do
    [[ -e $f ]] || die "missing expected file: $f"
done

# ---------------------------------------------------------------- transact workspace
if [[ ${GDCQC_SMOKE_KEEP:-0} == 1 && -d $WORK_BASE/OUT ]]; then
    log "KEEP=1: reusing existing workspace $WORK_BASE (skipping toy generation)"
else
    rm -rf "$WORK_BASE"
    mkdir -p "$OUT_DIR" "$SIM_INPUTS" "$RESOURCES_DIR"
    log "workspace prepared: $WORK_BASE"
fi

# ---------------------------------------------------------------- preflight tools
log "preflighting tools ..."

if command -v snakemake >/dev/null 2>&1; then
    SNAKEMAKE_BIN="$(command -v snakemake)"
else
    export PATH="/projects/standard/gdc/public/envs/snakemake/bin:$PATH"
    command -v snakemake >/dev/null 2>&1 || die "snakemake not found (tried envs/snakemake)"
    SNAKEMAKE_BIN="$(command -v snakemake)"
fi
log "  snakemake       : $SNAKEMAKE_BIN ($(snakemake --version))"

R_CANDIDATES=()
if command -v Rscript >/dev/null 2>&1; then
    R_CANDIDATES+=("$(command -v Rscript)")
fi
R_CANDIDATES+=(/common/software/install/manual/R/4.4.0-openblas-rocky8-fix/bin/Rscript)
RSCRIPT_BIN=""
for cand in "${R_CANDIDATES[@]}"; do
    if [[ -n $cand && -x $cand ]] && "$cand" -e 'suppressMessages(library(argparse)); suppressMessages(library(bigsnpr)); cat("ok\n")' >/dev/null 2>&1; then
        RSCRIPT_BIN="$cand"
        break
    fi
done
[[ -n $RSCRIPT_BIN ]] || die "no Rscript with argparse+bigsnpr found (tried: ${R_CANDIDATES[*]})"
export PATH="$(dirname "$RSCRIPT_BIN"):$PATH"
log "  Rscript         : $RSCRIPT_BIN"
"$RSCRIPT_BIN" --version | head -1 | sed 's/^/  /'

PLINK_BIN=""; PLINK2_BIN=""
for c in plink plink2; do
    if command -v $c >/dev/null 2>&1; then
        [[ $c == plink  ]] && PLINK_BIN="$(command -v $c)"
        [[ $c == plink2 ]] && PLINK2_BIN="$(command -v $c)"
    fi
done
[[ -n $PLINK_BIN  ]] || PLINK_BIN="/projects/standard/gdc/public/envs/plink/bin/plink"
[[ -n $PLINK2_BIN ]] || PLINK2_BIN="/projects/standard/gdc/public/envs/plink/bin/plink2"
[[ -x $PLINK_BIN && -x $PLINK2_BIN ]] || die "plink/plink2 binaries not found"
log "  plink           : $PLINK_BIN"
log "  plink2          : $PLINK2_BIN"

PRSCS_PYTHON=""
for cand in /projects/standard/gdc/public/envs/gdcPipeline/bin/python python3; do
    if [[ -x $cand ]] || command -v "$cand" >/dev/null 2>&1; then
        if "$cand" -c 'import numpy, scipy, h5py' >/dev/null 2>&1; then
            PRSCS_PYTHON="$(command -v "$cand" || echo "$cand")"
            break
        fi
    fi
done
[[ -n $PRSCS_PYTHON ]] || die "no python with numpy+scipy+h5py for PRS-CS (set PATH or export GDCQC env)"
log "  PRS-CS python   : $PRSCS_PYTHON"
"$PRSCS_PYTHON" -c 'import numpy, scipy, h5py; print("  (numpy %s scipy %s h5py %s)" % (numpy.__version__, scipy.__version__, h5py.__version__))'

# external PRS-CS resources (PRScsx.py and 1KG LD reference)
PRSCSX_CODE="/projects/standard/gdc/public/prs_methods/scripts/PRScsx"
PRSCSX_REF="/projects/standard/gdc/public/prs_methods/ref/ref_PRScsx/1kg_ref"
[[ -f $PRSCSX_CODE/PRScsx.py ]] || die "missing PRScsx.py at $PRSCSX_CODE"
[[ -d $PRSCSX_REF/ldblk_1kg_eur ]] || die "missing 1KG LD reference at $PRSCSX_REF"
log "  PRScsx.py code  : $PRSCSX_CODE"
log "  1KG LD ref      : $PRSCSX_REF"

# ---------------------------------------------------------------- toy generation
if [[ ! -e $METHOD_RUN/single_ct.done && ${GDCQC_SMOKE_KEEP:-0} != 1 ]]; then
    section "Fabricating real-SNP toy PLINK1 sets"
    runas "AFR toy" "$PRSCS_PYTHON" "$GENERATOR" "$SIM_INPUTS/$ANC1/study" "$NSNPS" "$NSAMPLES" "$SEED"
    runas "EUR toy" "$PRSCS_PYTHON" "$GENERATOR" "$SIM_INPUTS/$ANC2/study" "$NSNPS" "$NSAMPLES" "$((SEED+1))"
    for anc in $ANC1 $ANC2; do
        b="$SIM_INPUTS/$anc/study"
        log "  $anc toy: $(wc -l < "$b.bim" | tr -d ' ') snps, $(wc -l < "$b.fam" | tr -d ' ') samples @ $b.{bed,bim,fam}"
        head -3 "$b.bim"
    done
    # Held-out test sample (same variant set, different samples + unique "t"
    # IIDs, seed SEED+2) for the scoreTestPRS / score_test.sh evaluation
    # (Phase 3) so train/test are unambiguous.
    runas "test toy" "$PRSCS_PYTHON" "$GENERATOR" "$SIM_INPUTS/test/study" "$NSNPS" "$NSAMPLES" "$((SEED+2))" "t"
    b="$SIM_INPUTS/test/study"
    log "  test toy: $(wc -l < "$b.bim" | tr -d ' ') snps, $(wc -l < "$b.fam" | tr -d ' ') samples @ $b.{bed,bim,fam}"
else
    log "workspace reused, skipping toy generation"
fi

# ---------------------------------------------------------------- seed the sim subgraph
# simulateBivariatePhenotypes' bfile branch hardcodes process-local CTSLEB scratch
# paths, and the QC->initialFilter.pgen graph is ambiguous -- bypass both by
# pre-creating its declared outputs; they must be NEWER than the toy inputs.
section "Seeding phenotypeSimulation outputs (bypass sim + QC subgraphs)"
SIM_OUT="$OUT_DIR/simulations/${ANC1}_${ANC2}"
mkdir -p "$SIM_OUT"
for anc in $ANC1 $ANC2; do
    sleep 1
    for ext in bed bim fam; do
        cp "$SIM_INPUTS/$anc/study.$ext" "$SIM_OUT/${anc}_simulation.$ext"
    done
done
touch "$SIM_OUT"
ls -la "$SIM_OUT"

# ---------------------------------------------------------------- write config
section "Writing test config -> $CONFIG"
cat > "$CONFIG" <<EOF
OUT_DIR: $(printf %q "$OUT_DIR")
INPUT: /projects/standard/gdc/public/toyData/allSubjects_hg38_chr{CHR}.vcf.gz
REF: /projects/standard/gdc/public/Ref
chromosomes:
  - 20
  - 21

phenotypeSimulation:
  ancestries: ["$ANC1", "$ANC2"]
  input_prefixes:
    ${ANC1}: $(printf %q "$SIM_INPUTS/$ANC1/study")
    ${ANC2}: $(printf %q "$SIM_INPUTS/$ANC2/study")

prsPipeline:
  path_plink: $(printf %q "$PLINK_BIN")
  path_plink2: $(printf %q "$PLINK2_BIN")
  n_total_gwas: ${NSAMPLES}
  binary_target: "F"
  seed: ${SEED}
  gwas_fraction: 0.5
  phenotype_index: 1
  test_sample: $(printf %q "$SIM_INPUTS/test/study")

prsMethods:
  resource_dir: $(printf %q "$RESOURCES_DIR")
  containers:
    prsv2: null
    singleprshelper: null
  single_prscs:
    path_code: $(printf %q "$PRSCSX_CODE")
    ld_ref_dir: $(printf %q "$PRSCSX_REF")
    seed: ${SEED}
    path_python: $(printf %q "$PRSCS_PYTHON")
  multi_prscsx:
    seed: ${SEED}
EOF
sed -n '1,2p;s/^\(IN\|OUT\|prs\|INPUT\|REF\)/  \1/p' "$CONFIG" 2>/dev/null | head -5
log "config written ($(wc -l < "$CONFIG") lines)"

# ---------------------------------------------------------------- dry run
section "Dry-run of run_singleAncestryPRSPipelines + runMultiAncestryPRSCSx"
(
    cd "$WORKFLOW"
    snakemake -n -j "$CPUS" $TRACE_ARGS run_singleAncestryPRSPipelines --configfile "$CONFIG" 2>&1 \
        | tee "$WORK_BASE/dryrun.log" \
        | grep -E "^(Job stats|runSingleAncestry|total|waiting|Nothing to be done|Job counts)" \
        || true
)
(
    cd "$WORKFLOW"
    snakemake -n -j "$CPUS" $TRACE_ARGS runMultiAncestryPRSCSx --configfile "$CONFIG" 2>&1 \
        | tee "$WORK_BASE/dryrun_mprscsx.log" \
        | grep -E "^(Job stats|runMultiAncestry|prepare|total|waiting|Nothing to be done|Job counts)" \
        || true
)
if [[ ${GDCQC_SMOKE_DRYRUN:-0} == 1 ]]; then
    log "GDCQC_SMOKE_DRYRUN=1, stopping after dry-run"
    exit 0
fi

# ---------------------------------------------------------------- real run
section "Running the single-ancestry PRS pipelines (snakemake -j $CPUS)"
if [[ -e $METHOD_RUN/single_ct.done && -e $METHOD_RUN/single_prscs.done && ${GDCQC_SMOKE_KEEP:-0} == 1 ]]; then
    log "KEEP=1 and all done markers present -- skipping execution"
else
    (
        cd "$WORKFLOW"
        snakemake -j "$CPUS" $TRACE_ARGS \
            --verbose \
            run_singleAncestryPRSPipelines \
            --configfile "$CONFIG" \
            > "$RUN_LOG" 2>&1
    )
    SNAKE_RC=$?
    log "snakemake finished with exit code $SNAKE_RC (see $RUN_LOG)"
    if [[ $SNAKE_RC -ne 0 ]]; then
        log "--- last 60 lines of $RUN_LOG ---"
        tail -60 "$RUN_LOG"
        die "snakemake run failed (rc=$SNAKE_RC)"
    fi
fi

# ---------------------------------------------------------------- test evaluation (Phase 3)
section "Running scoreTestPRS on the held-out test sample (snakemake -j $CPUS)"
SCORE_LOG="$WORK_BASE/scoring_run.log"
if [[ -e $METHOD_RUN/scoreTestPRS.done && ${GDCQC_SMOKE_KEEP:-0} == 1 ]]; then
    log "KEEP=1 and scoreTestPRS.done present -- skipping execution"
else
    (
        cd "$WORKFLOW"
        snakemake -j "$CPUS" $TRACE_ARGS \
            --verbose \
            run_scoreTestPRS \
            --configfile "$CONFIG" \
            > "$SCORE_LOG" 2>&1
    )
    SNAKE_RC=$?
    log "score test finished with exit code $SNAKE_RC (see $SCORE_LOG)"
    if [[ $SNAKE_RC -ne 0 ]]; then
        log "--- last 60 lines of $SCORE_LOG ---"
        tail -60 "$SCORE_LOG"
        die "score test run failed (rc=$SNAKE_RC)"
    fi
fi

# ---------------------------------------------------------------- multi-ancestry PRS-CSx
section "Running joint multi-ancestry PRS-CSx (snakemake -j $CPUS)"
MPRSCSX_LOG="$WORK_BASE/mprscsx_run.log"
if [[ -e $METHOD_RUN/multi_prscsx.done && ${GDCQC_SMOKE_KEEP:-0} == 1 ]]; then
    log "KEEP=1 and multi_prscsx.done present -- skipping execution"
else
    (
        cd "$WORKFLOW"
        snakemake -j "$CPUS" $TRACE_ARGS \
            --verbose \
            runMultiAncestryPRSCSx \
            --configfile "$CONFIG" \
            > "$MPRSCSX_LOG" 2>&1
    )
    SNAKE_RC=$?
    log "multi_prscsx run finished with exit code $SNAKE_RC (see $MPRSCSX_LOG)"
    if [[ $SNAKE_RC -ne 0 ]]; then
        log "--- last 60 lines of $MPRSCSX_LOG ---"
        tail -60 "$MPRSCSX_LOG"
        die "multi_prscsx run failed (rc=$SNAKE_RC)"
    fi
fi

# ---------------------------------------------------------------- verification
section "Verifying method outputs"
FAIL=0
check() {  # check <label> <relative-path-to-method-run> [min-expected-size]
    local label="$1" rel="$2" min="$3"
    local path="$METHOD_RUN/$rel"
    if [[ -s $path ]] && [[ $(stat -c %s "$path") -ge ${min:-1} ]]; then
        printf '  PASS  %-12s %s (%d bytes)\n' "$label" "$rel" "$(stat -c %s "$path")"
    else
        printf '  FAIL  %-12s %s MISSING/EMPTY\n' "$label" "$rel"
        FAIL=$((FAIL+1))
    fi
}
donef() {  # donef <method>
    if [[ -e $METHOD_RUN/$1.done ]]; then
        printf '  PASS  %-12s %s.done\n' "$1" "$1"
    else
        printf '  FAIL  %-12s %s.done MISSING\n' "$1" "$1"
        FAIL=$((FAIL+1))
    fi
}

for m in single_ct single_prsice single_prscs single_ldpred2 single_lassosum2 multi_prscsx; do
    donef "$m"
done

check "CT"        "single_ct/CT/CT_prs_results.txt"            100
check "PRSice2"   "single_prsice/PRSice2/prs_method/PRSice2_outputs.prsice" 100
check "PRS-CS"    "single_prscs/prs_pipeline/PRScs/PRScs_${ANC1}_combined_weights.txt" 100
check "PRS-CS"    "single_prscs/prs_pipeline/PRScs/PRScs_${ANC1}_score.sscore" 500
check "PRS-CS"    "single_prscs/prs_pipeline/PRScs/${ANC1}_PRS_sscore_Rsqr.txt" 1
check "LDpred2"   "single_ldpred2/prs_method_individual_scores.txt" 500
check "LDpred2"   "single_ldpred2/prs_method_performance.csv" 10
check "lassosum2" "single_lassosum2/prs_method_final_res.txt" 10
check "lassosum2" "single_lassosum2/prs_method_full_predictions.csv" 500

# --- joint multi-ancestry PRS-CSx ---
check "mPRS-CSx"  "multi_prscsx/prs_pipeline/PRScsx/PRScsx_${ANC1}_combined_weights.txt" 100
check "mPRS-CSx"  "multi_prscsx/prs_pipeline/PRScsx/PRScsx_${ANC2}_combined_weights.txt" 100
check "mPRS-CSx"  "multi_prscsx/prs_pipeline/PRScsx/PRScsx_joint_${ANC1}_score.sscore" 500
check "mPRS-CSx"  "multi_prscsx/prs_pipeline/PRScsx/PRScsx_joint_${ANC2}_score.sscore" 500
check "mPRS-CSx"  "multi_prscsx/prs_pipeline/PRScsx/${ANC1}_PRS_sscore_Rsqr.txt" 1
check "mPRS-CSx"  "multi_prscsx/prs_pipeline/PRScsx/${ANC2}_PRS_sscore_Rsqr.txt" 1

# --- Phase 3: held-out test evaluation (score_test.sh) ---
for m in CT LDpred2_inf LDpred2_grid lassosum2 PRSice2 PRScsx_${ANC1}; do
    check "test-eval" "test_evaluation/${m}_results.txt" 10
    check "test-eval" "test_evaluation/${m}_scores.txt" 10
done
donef "scoreTestPRS"

echo
log "Key result snippets:"
[[ -f $METHOD_RUN/single_prscs/prs_pipeline/PRScs/${ANC1}_PRS_sscore_Rsqr.txt ]] \
    && { echo "  PRS-CS Rsqr:"; sed 's/^/    /' "$METHOD_RUN/single_prscs/prs_pipeline/PRScs/${ANC1}_PRS_sscore_Rsqr.txt"; }
[[ -s $METHOD_RUN/single_prsice/PRSice2/prs_method/PRSice2_outputs.prsice ]] \
    && { echo "  PRSice top row:"; head -2 "$METHOD_RUN/single_prsice/PRSice2/prs_method/PRSice2_outputs.prsice" | sed 's/^/    /'; }
[[ -s $METHOD_RUN/single_ldpred2/prs_method_performance.csv ]] \
    && { echo "  LDpred2 performance:"; sed 's/^/    /' "$METHOD_RUN/single_ldpred2/prs_method_performance.csv"; }
echo "  Multi-prscsx R2 (ANC1=${ANC1} / ANC2=${ANC2}):"
for anc in $ANC1 $ANC2; do
    f="$METHOD_RUN/multi_prscsx/prs_pipeline/PRScsx/${anc}_PRS_sscore_Rsqr.txt"
    [[ -s $f ]] && { echo -n "    $anc: "; sed 's/^/    /' "$f"; }
done
echo "  Held-out test R2 (score_test.sh, test_evaluation/):"
for m in CT LDpred2_inf LDpred2_grid lassosum2 PRSice2 PRScsx_${ANC1}; do
    f="$METHOD_RUN/test_evaluation/${m}_results.txt"
    [[ -s $f ]] && { echo -n "    $m: "; sed 's/^/    /' "$f" | grep -Ei "r2|r.sq|p.value" | head -1; }
done

section "SUMMARY"
ELAPSED="$(( $(date +%s) - START_EPOCH ))"
if [[ $FAIL -eq 0 ]]; then
    echo "  ALL SINGLE-ANCESTRY PRS METHODS + MULTI_PRSCSX + HELD-OUT TEST EVALUATION PASSED in ${ELAPSED}s"
    echo "  outputs under: $METHOD_RUN"
else
    echo "  $FAIL VERIFICATION CHECK(S) FAILED (elapsed ${ELAPSED}s)"
    echo "  see:  $RUN_LOG (methods) / $SCORE_LOG (score test) / $MPRSCSX_LOG (multi prscsx)"
    echo "  logs: $(cd "$OUT_DIR/logs" && ls -1 2>/dev/null | tr '\n' ' ')"
    exit 1
fi

echo
echo "  Directories of interest:"
echo "    toy inputs : $SIM_INPUTS"
echo "    resources  : $RESOURCES_DIR"
echo "    outputs    : $PRS_OUT"
echo "    run log    : $RUN_LOG"
echo "    score log  : $SCORE_LOG"
echo "    mprscsx log: $MPRSCSX_LOG"
echo
echo "  Job stats from dry run:"
grep -E "^(runSingleAncestry|total|waiting)" "$WORK_BASE/dryrun.log" || true
echo
echo "  Next up: remaining multi-ancestry methods (multi_ctsleb/multi_ldpred2/multi_sdprs/multi_prosper) + PROSPER ref_bim.txt blocker."
exit 0