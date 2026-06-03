#!/bin/bash
# Single-Ancestry PRS Pipeline
#
# Runs single-ancestry PRS methods (PRSice2, LDPred2) on ancestry-stratified data.
#
# Usage: bash run_prs.sh [msi|sandbox|other|local]
#
# Prerequisites:
# - Run setup_env.sh first (or source it)
# - Ancestry labels available
# - Summary statistics and phenotype data prepared

ENV="${1:-msi}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$SCRIPT_DIR/setup_env.sh" "$ENV"

# Step 1: Create configuration file
mkdir -p ~/prs_lab
cd ~/prs_lab

cat > config_prs.yaml << 'EOF'
INPUT: "/path/to/data/chr{CHR}.vcf.gz"
REF: "/path/to/reference/storage"
OUT_DIR: "/path/to/output/directory"
local-storage-prefix: "/path/to/.snakemake/storage"

chromosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]

ancestry:
    ancestry_file: "/path/to/ancestry_labels.tsv"

prsMethods:
    resource_dir: "/path/to/prs_resources"
    single_prsice:
      enabled: true
    single_ldpred2:
      enabled: true
    single_ct:
      enabled: false
    single_prscs:
      enabled: false
    single_lassosum2:
      enabled: false

conda-frontend: mamba
EOF

echo "Config written to config_prs.yaml"
echo "Edit paths in config_prs.yaml before running."

# Step 2: Run all enabled PRS methods
cd GDCGenomicsQC/workflow
echo "Step 2: Running all enabled single-ancestry PRS methods..."
gdcgenomicsqc --configfile ../config_prs.yaml runAllEnabledPRS -j 4

# Optional: Run individual methods
echo ""
echo "To run individual methods:"
echo "  gdcgenomicsqc --configfile ../config_prs.yaml runSingleAncestryPRSice -j 4"
echo "  gdcgenomicsqc --configfile ../config_prs.yaml runSingleAncestryCT -j 4"

echo "PRS pipeline complete."
echo "Outputs in: PRS_OUT_DIR/method_runs/"
