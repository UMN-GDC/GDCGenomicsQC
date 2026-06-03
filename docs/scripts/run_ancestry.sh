#!/bin/bash
# Global Ancestry Classification Pipeline
#
# Runs ancestry classification using PCA (or UMAP) with Random Forest
# on reference-projected coordinates. Produces posterior probabilities,
# ancestry classifications, and confusion matrices.
#
# Usage: bash run_ancestry.sh [msi|sandbox|other|local]
#
# Prerequisites:
# - Run setup_env.sh first (or source it)
# - QC-filtered genotype data available

ENV="${1:-msi}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$SCRIPT_DIR/setup_env.sh" "$ENV"

# Step 1: Create configuration file
mkdir -p ~/ancestry_lab
cd ~/ancestry_lab

cat > config_ancestry.yaml << 'EOF'
INPUT: "/path/to/data/chr{CHR}.vcf.gz"
OUT_DIR: "/path/to/output/directory"
REF: "/path/to/reference/data"
local-storage-prefix: "/path/to/.snakemake/storage"

chromosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]

ancestry:
    threshold: 0.8
    model: "pca"

relatedness:
    method: "king"
    king_cutoff: 0.0884

internalPCA:
    method: "plink2"
    npc: 20

localAncestry:
    RFMIX: true
    test: true
    thin_subjects: 0.1
    figures: "figures"
    chromosomes: null

thin: false
conda-frontend: mamba
EOF

echo "Config written to config_ancestry.yaml"
echo "Edit INPUT, OUT_DIR, REF in config_ancestry.yaml before running."

# Step 2: Run classification pipeline
cd GDCGenomicsQC/workflow
echo "Step 2: Running classification pipeline..."
gdcgenomicsqc --configfile ../config_ancestry.yaml classifyAncestry -j 10

# Step 3: Compare models (PCA vs UMAP)
echo "Step 3: Comparing PCA vs UMAP models..."
gdcgenomicsqc --configfile ../config_ancestry.yaml classifyAncestry

# Step 4: Ancestry-specific subsetting
echo "Step 4: Extracting ancestry-specific data..."
gdcgenomicsqc --configfile ../config_ancestry.yaml convertNfilt/CHR=20/subset=EUR

echo "Ancestry classification complete."
echo "Outputs in: OUT_DIR/01-globalAncestry/"
