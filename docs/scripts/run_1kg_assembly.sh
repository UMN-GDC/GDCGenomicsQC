#!/bin/bash
# 1000 Genomes Reference Assembly Pipeline
#
# Downloads, processes, and assembles the 1000 Genomes high-coverage
# reference panel into PLINK2 format for ancestry classification and QC.
#
# Usage: bash run_1kg_assembly.sh [msi|sandbox|other|local]
#
# Prerequisites:
# - Run setup_env.sh first (or source it)
# - ~50GB available for the reference panel

ENV="${1:-msi}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$SCRIPT_DIR/setup_env.sh" "$ENV"

# Step 1: Configure reference paths
mkdir -p ~/reference_lab
cd ~/reference_lab

cat > config_reference.yaml << 'EOF'
REF: "/path/to/reference/storage"
OUT_DIR: "/path/to/output/directory"
local-storage-prefix: "/path/to/.snakemake/storage"

chromosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]

conda-frontend: mamba
EOF

echo "Config written to config_reference.yaml"
echo "Edit REF, OUT_DIR, and local-storage-prefix in config_reference.yaml before running."

# Step 2: Download metadata
echo "Step 2: Downloading metadata..."
cd GDCGenomicsQC/workflow
gdcgenomicsqc --configfile ../config_reference.yaml kgMeta -j 4

# Step 3: Download VCF data (22 chromosomes in parallel)
echo "Step 3: Downloading VCF data..."
gdcgenomicsqc --configfile ../config_reference.yaml kgData -j 22

# Step 4: Assemble into PLINK2 format
echo "Step 4: Assembling into PLINK2 format..."
gdcgenomicsqc --configfile ../config_reference.yaml kgAssemble -j 8

echo "Reference assembly complete."
