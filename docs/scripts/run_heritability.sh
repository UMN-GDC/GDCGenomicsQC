#!/bin/bash
# SNP Heritability Estimation Pipeline
#
# Estimates SNP heritability using various methods (AdjHE, GCTA, PredLMM, SWD)
# via the MASH framework.
#
# Usage: bash run_heritability.sh [msi|sandbox|other|local]
#
# Prerequisites:
# - Run setup_env.sh first (or source it)
# - Ancestry labels available (or ancestry classification completed)
# - Phenotype and covariate files prepared

ENV="${1:-msi}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$SCRIPT_DIR/setup_env.sh" "$ENV"

# Step 1: Create configuration file
mkdir -p ~/heritability_lab
cd ~/heritability_lab

cat > config_heritability.yaml << 'EOF'
INPUT: "/path/to/data/chr{CHR}.vcf.gz"
REF: "/path/to/reference/storage"
OUT_DIR: "/path/to/output/directory"
local-storage-prefix: "/path/to/.snakemake/storage"

chromosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]

ancestry:
    ancestry_file: "/path/to/ancestry_labels.tsv"

snpHerit:
    pheno: "/path/to/phenotype.tsv"
    covar: "/path/to/covariates.tsv"
    method: "AdjHE"
    npc: 10
    mpheno: "BMI"
    loop_covars: false
    Naive: false

conda-frontend: mamba
EOF

echo "Config written to config_heritability.yaml"
echo "Edit INPUT, REF, OUT_DIR, ancestry, snpHerit paths before running."

# Step 2: Run heritability estimation
cd GDCGenomicsQC/workflow
echo "Step 2: Running SNP heritability estimation..."
gdcgenomicsqc --configfile ../config_heritability.yaml snpHerit -j 4

echo "Heritability estimation complete."
echo "Outputs in: OUT_DIR/03-snpHeritability/"
