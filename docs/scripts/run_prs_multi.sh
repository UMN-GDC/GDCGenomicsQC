#!/bin/bash
# Multi-Ancestry PRS Pipeline
#
# Runs multi-ancestry PRS methods (CT-SLeB Multi, PRScsx, SDPRSx)
# using summary statistics from multiple ancestries.
#
# Usage: bash run_prs_multi.sh [msi|sandbox|other|local]
#
# Prerequisites:
# - Run setup_env.sh first (or source it)
# - Ancestry classification completed
# - Training and target summary statistics prepared
# - LD reference panels available for PRScsx and SDPRSx

ENV="${1:-msi}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$SCRIPT_DIR/setup_env.sh" "$ENV"

# Step 1: Create configuration file
mkdir -p ~/prs_multi_lab
cd ~/prs_multi_lab

cat > config_prs_multi.yaml << 'EOF'
INPUT: "/path/to/data/chr{CHR}.vcf.gz"
REF: "/path/to/reference/storage"
OUT_DIR: "/path/to/output/directory"
PRS_OUT_DIR: "/path/to/prs/output"
local-storage-prefix: "/path/to/.snakemake/storage"

chromosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]

ancestry:
    model: "pca"
    threshold: 0.8

prsMethods:
    resource_dir: "/path/to/prs_resources"
    multi_ctsleb:
      enabled: true
    multi_prscsx:
      enabled: true
      ld_ref_dir: "/path/to/prscsx/ld/ref"
    multi_sdprs:
      enabled: true
      ld_ref_dir: "/path/to/sdprs/ld"
    multi_ldpred2:
      enabled: false
    multi_prosper:
      enabled: false
    single_ct:
      enabled: false
    single_prsice:
      enabled: false
    single_prscs:
      enabled: false
    single_ldpred2:
      enabled: false
    single_lassosum2:
      enabled: false

conda-frontend: mamba
EOF

echo "Config written to config_prs_multi.yaml"
echo "Edit paths in config_prs_multi.yaml before running."

# Step 2: Run all enabled multi-ancestry PRS methods
cd GDCGenomicsQC/workflow
echo "Step 2: Running all enabled multi-ancestry PRS methods..."
gdcgenomicsqc --configfile ../config_prs_multi.yaml runAllEnabledPRS -j 4

# Optional: Run individual methods
echo ""
echo "To run individual methods:"
echo "  gdcgenomicsqc --configfile ../config_prs_multi.yaml runMultiAncestryPRSCSx -j 4"
echo "  gdcgenomicsqc --configfile ../config_prs_multi.yaml runMultiAncestrySDPRS -j 4"
echo "  gdcgenomicsqc --configfile ../config_prs_multi.yaml runMultiAncestryCTSLEB -j 4"

echo "Multi-ancestry PRS pipeline complete."
echo "Outputs in: PRS_OUT_DIR/method_runs/"
