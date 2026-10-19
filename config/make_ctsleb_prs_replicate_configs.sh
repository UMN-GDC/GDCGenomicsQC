#!/usr/bin/env bash
# Create one PRS configuration per simulated phenotype without re-running simulation.
set -euo pipefail

CONFIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEMPLATE="${1:-${CONFIG_DIR}/msiCTSLEBPRS_corrected.yaml}"
N_REPLICATES="${2:-10}"

[[ -f "$TEMPLATE" ]] || { echo "Missing template: $TEMPLATE" >&2; exit 1; }
[[ "$N_REPLICATES" =~ ^[1-9][0-9]*$ ]] || { echo "N_REPLICATES must be a positive integer" >&2; exit 1; }
OUT_DIR=$(awk -F'"' '/^OUT_DIR:/ {print $2; exit}' "$TEMPLATE")
[[ -n "$OUT_DIR" ]] || { echo "Could not read OUT_DIR from $TEMPLATE" >&2; exit 1; }
PRS_ROOT="${OUT_DIR}/prs_inputs_ctsleb/AFR_EUR"
CONFIG_STEM="$(basename "$TEMPLATE" .yaml)"

for ((replicate = 1; replicate <= N_REPLICATES; replicate++)); do
    tag=$(printf "rep%02d" "$replicate")
    output="${CONFIG_DIR}/${CONFIG_STEM}_${tag}.yaml"

    sed \
        -e "s|${PRS_ROOT}|${PRS_ROOT}/${tag}|g" \
        -e "s|phenotype_index: 1|phenotype_index: ${replicate}|" \
        "$TEMPLATE" > "$output"

    echo "Wrote ${output} for phenotype ${replicate}"
done
