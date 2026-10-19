#!/usr/bin/env bash
# Run the CTSLEB phenotype simulation and PRS benchmark on an Agate login node.
#
# The script is resumable: completed Snakemake outputs are reused, while failed
# or partial jobs are retried with --rerun-incomplete.  Do not run two copies
# of this script against the same corrected output directory at the same time.

set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-/scratch.global/saonli/GDCGenomicsQC}"
WORKFLOW_DIR="${WORKFLOW_DIR:-${PROJECT_DIR}/workflow}"
CONFIG_DIR="${CONFIG_DIR:-${PROJECT_DIR}/config}"
BASE_CONFIG="${BASE_CONFIG:-${CONFIG_DIR}/msiCTSLEBPRS_manuscript.yaml}"
PROFILE_DIR="${PROFILE_DIR:-${PROJECT_DIR}/profiles/hpc}"
N_REPLICATES="${N_REPLICATES:-10}"

# Methods enabled by default. CT-SLEB uses the corrected simulated BIM files,
# which passed the unique-variant-ID preflight check for AFR and EUR.
RUN_SINGLE="${RUN_SINGLE:-1}"
RUN_PRSCSX="${RUN_PRSCSX:-1}"
RUN_PROSPER="${RUN_PROSPER:-1}"
RUN_CTSLEB="${RUN_CTSLEB:-1}"

# PRS-CSx needs more than the rule's historical 12-hour/32-GB default on the
# full CTSLEB data.  Values are Slurm minutes and MB.
PRSCX_MEM_MB="${PRSCX_MEM_MB:-64000}"
PRSCX_RUNTIME_MIN="${PRSCX_RUNTIME_MIN:-2160}"
PROSPER_MEM_MB="${PROSPER_MEM_MB:-64000}"
PROSPER_RUNTIME_MIN="${PROSPER_RUNTIME_MIN:-2160}"
CTSLEB_MEM_MB="${CTSLEB_MEM_MB:-64000}"
CTSLEB_RUNTIME_MIN="${CTSLEB_RUNTIME_MIN:-1440}"

usage() {
    cat <<'EOF'
Usage: bash run_ctsleb_prs_benchmark.sh

Run this script from an Agate login node. Configure it with environment
variables instead of editing it, for example:

  N_REPLICATES=1 bash run_ctsleb_prs_benchmark.sh
  N_REPLICATES=10 bash run_ctsleb_prs_benchmark.sh
  RUN_CTSLEB=0 N_REPLICATES=1 bash run_ctsleb_prs_benchmark.sh

Defaults run C+T, LDpred2, lassosum2, PRS-CSx, PROSPER, and CT-SLEB for all
ten simulated phenotypes. Set RUN_CTSLEB=0 to omit CT-SLEB. Results are collected in
output_ctsleb_prs_manuscript/benchmark_results_raw.tsv.
EOF
}

case "${1:-}" in
    -h|--help) usage; exit 0 ;;
    "") ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
esac

[[ "${N_REPLICATES}" =~ ^[1-9][0-9]*$ ]] || {
    echo "N_REPLICATES must be a positive integer." >&2
    exit 2
}
[[ -f "${BASE_CONFIG}" ]] || { echo "Missing base config: ${BASE_CONFIG}" >&2; exit 2; }
[[ -d "${WORKFLOW_DIR}" ]] || { echo "Missing workflow directory: ${WORKFLOW_DIR}" >&2; exit 2; }
[[ -d "${PROFILE_DIR}" ]] || { echo "Missing Snakemake profile: ${PROFILE_DIR}" >&2; exit 2; }

# The controller must run on a login node, not inside a Slurm allocation.
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    echo "Run this driver from an Agate login node, not within Slurm job ${SLURM_JOB_ID}." >&2
    exit 2
fi

# MSI's Conda activation hooks inspect optional variables that are unset in a
# strict Bash shell. Disable only the unset-variable guard during activation.
set +u
source /projects/standard/gdc/public/envs/load_miniconda3.sh
conda activate /projects/standard/gdc/public/envs/snakemake
set -u

OUT_DIR=$(awk -F'"' '/^OUT_DIR:/ {print $2; exit}' "${BASE_CONFIG}")
[[ -n "${OUT_DIR}" ]] || { echo "Could not read OUT_DIR from ${BASE_CONFIG}" >&2; exit 2; }
CONFIG_STEM="$(basename "${BASE_CONFIG}" .yaml)"
SIM_MANIFEST="${OUT_DIR}/simulations/AFR_EUR/simulation_truth/simulation_manifest.tsv"
RUN_MANIFEST="${OUT_DIR}/benchmark_run_manifest.tsv"
CONTROLLER_LOG_DIR="${OUT_DIR}/benchmark_controller_logs"
COMPUTE_JOB_IDS="${OUT_DIR}/benchmark_compute_jobids.tsv"

run_snakemake() {
    local config_file="$1"
    shift
    local stage="${BENCHMARK_STAGE:-controller}"
    local controller_log="${CONTROLLER_LOG_DIR}/${stage}_$(date +%Y%m%dT%H%M%S).log"
    local status
    mkdir -p "${CONTROLLER_LOG_DIR}"
    set +e
    snakemake -p \
        --profile "${PROFILE_DIR}" \
        --configfile "${config_file}" \
        --rerun-incomplete \
        --latency-wait 180 \
        "$@" 2>&1 | tee "${controller_log}"
    status=${PIPESTATUS[0]}
    set -e
    { grep -Eo "SLURM jobid [0-9]+" "${controller_log}" || true; } | awk '{print $3}' | sort -u | while read -r job_id; do
        [[ -n "${job_id}" ]] || continue
        printf '%s\t%s\t%s\n' "${stage}" "${job_id}" "${controller_log}" >> "${COMPUTE_JOB_IDS}"
    done
    return "${status}"
}

cd "${WORKFLOW_DIR}"
mkdir -p "${OUT_DIR}"
if [[ ! -s "${COMPUTE_JOB_IDS}" ]]; then
    printf 'benchmark_stage\tjob_id\tcontroller_log\n' > "${COMPUTE_JOB_IDS}"
fi

if [[ ! -s "${SIM_MANIFEST}" ]]; then
    echo "[$(date)] Simulating the phenotype replicates configured in ${BASE_CONFIG}."
    BENCHMARK_STAGE="simulation"
    run_snakemake "${BASE_CONFIG}" simulateBivariatePhenotypes
else
    echo "[$(date)] Reusing existing phenotype simulation: ${SIM_MANIFEST}"
fi

bash "${CONFIG_DIR}/make_ctsleb_prs_replicate_configs.sh" "${BASE_CONFIG}" "${N_REPLICATES}"

mkdir -p "${OUT_DIR}"
if [[ ! -s "${RUN_MANIFEST}" ]]; then
    printf 'timestamp\treplicate\tstage\tstatus\tconfig\n' > "${RUN_MANIFEST}"
fi

for ((replicate = 1; replicate <= N_REPLICATES; replicate++)); do
    tag=$(printf 'rep%02d' "${replicate}")
    config_file="${CONFIG_DIR}/${CONFIG_STEM}_${tag}.yaml"
    targets=()
    resources=()

    (( RUN_SINGLE == 1 )) && targets+=(run_singleAncestryPRS)
    if (( RUN_PRSCSX == 1 )); then
        targets+=(runMultiAncestryPRSCSx)
        resources+=(--set-resources "runMultiAncestryPRSCSx:mem_mb=${PRSCX_MEM_MB}")
        resources+=(--set-resources "runMultiAncestryPRSCSx:runtime=${PRSCX_RUNTIME_MIN}")
    fi
    if (( RUN_PROSPER == 1 )); then
        targets+=(runMultiAncestryPROSPER)
        resources+=(--set-resources "runMultiAncestryPROSPER:mem_mb=${PROSPER_MEM_MB}")
        resources+=(--set-resources "runMultiAncestryPROSPER:runtime=${PROSPER_RUNTIME_MIN}")
    fi
    if (( RUN_CTSLEB == 1 )); then
        targets+=(runMultiAncestryCTSLEB)
        resources+=(--set-resources "runMultiAncestryCTSLEB:mem_mb=${CTSLEB_MEM_MB}")
        resources+=(--set-resources "runMultiAncestryCTSLEB:runtime=${CTSLEB_RUNTIME_MIN}")
    fi

    if (( ${#targets[@]} == 0 )); then
        echo "No methods selected. Set at least one RUN_* variable to 1." >&2
        exit 2
    fi

    BENCHMARK_STAGE="${tag}"
    printf '%s\t%s\t%s\t%s\t%s\n' "$(date --iso-8601=seconds)" "${tag}" \
        "submit" "started" "${config_file}" >> "${RUN_MANIFEST}"
    echo "[$(date)] Starting ${tag}: ${targets[*]}"

    # Parse resource overrides before the end-of-options marker; target rules
    # after -- cannot be consumed as resource assignments.
    if run_snakemake "${config_file}" "${resources[@]}" -- "${targets[@]}"; then
        printf '%s\t%s\t%s\t%s\t%s\n' "$(date --iso-8601=seconds)" "${tag}" \
            "submit" "completed" "${config_file}" >> "${RUN_MANIFEST}"
    else
        printf '%s\t%s\t%s\t%s\t%s\n' "$(date --iso-8601=seconds)" "${tag}" \
            "submit" "failed" "${config_file}" >> "${RUN_MANIFEST}"
        echo "${tag} failed. Fix the reported rule, then rerun this same command to resume." >&2
        exit 1
    fi
done

python3 "${CONFIG_DIR}/collect_ctsleb_prs_results.py" \
    --output-dir "${OUT_DIR}" \
    --replicates "${N_REPLICATES}"
python3 "${CONFIG_DIR}/collect_prs_compute_metrics.py" \
    --workflow-dir "${WORKFLOW_DIR}" \
    --output-dir "${OUT_DIR}" \
    --job-id-file "${COMPUTE_JOB_IDS}"

echo "[$(date)] Benchmark complete."
echo "Raw result inventory: ${OUT_DIR}/benchmark_results_raw.tsv"
echo "Compute-time inventory: ${OUT_DIR}/benchmark_compute_metrics.tsv"
