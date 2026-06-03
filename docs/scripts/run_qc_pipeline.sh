#!/bin/bash
# Basic QC Pipeline
#
# Runs the two-stage QC pipeline: Initial QC (missingness, LD pruning)
# and Standard QC (MAF, HWE, heterozygosity, sex check).
#
# Usage: bash run_qc_pipeline.sh [msi|sandbox|other|local]
#
# Prerequisites:
# - Run setup_env.sh first (or source it)
# - Input genotype data available at the path specified in config

ENV="${1:-msi}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$SCRIPT_DIR/setup_env.sh" "$ENV"

# Step 1: Create configuration file
mkdir -p ~/qc_lab
cd ~/qc_lab

cat > config_qc.yaml << 'EOF'
INPUT: "/path/to/data/chr{CHR}.vcf.gz"
OUT_DIR: "/path/to/output/directory"
REF: "/path/to/reference/data"
local-storage-prefix: "/path/to/.snakemake/storage"

chromosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]

relatedness:
    method: "king"
    king_cutoff: 0.0884

SEX_CHECK: true
thin: false
conda-frontend: mamba

internalPCA:
    method: "plink2"
    npc: 20
EOF

echo "Config written to config_qc.yaml"
echo "Edit INPUT, OUT_DIR, REF in config_qc.yaml before running."

# Step 2: Run Initial QC
cd GDCGenomicsQC/workflow
echo "Step 2: Running Initial QC..."
gdcgenomicsqc --configfile ../config_qc.yaml full/f1.pgen -j 10

# Step 3: Run Standard QC
echo "Step 3: Running Standard QC..."
gdcgenomicsqc --configfile ../config_qc.yaml full/f1.b38.f2.pgen -j 10

# Step 4: Run ancestry-specific QC (if ancestry labels available)
echo "Step 4: Running ancestry-specific QC..."
gdcgenomicsqc --configfile ../config_qc.yaml EUR/f1.b38.f2.pgen -j 10

echo "QC pipeline complete."
