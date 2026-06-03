#!/bin/bash
# Setup environment for GDCGenomicsQC pipeline
#
# Usage: source setup_env.sh [msi|sandbox|other|local]
#
# Select your environment:
#   msi     - MSI HPC Agate
#   sandbox - Sandbox environment
#   other   - Other HPC (customize module path)
#   local   - Local Snakemake (conda only)

ENV="${1:-msi}"

case "$ENV" in
  msi)
    echo "Setting up MSI HPC environment..."
    module use /projects/standard/gdc/public/GDCGenomicsQC/envs
    module load gdcgenomicsMSI
    conda activate snakemake
    ;;
  sandbox)
    echo "Setting up Sandbox environment..."
    module use /scratch.global/GDC/GDCGenomicsQC/envs
    module load gdcgenomicsSandbox
    conda activate snakemake
    ;;
  other)
    echo "Setting up Other HPC environment..."
    echo "Edit MODULE_PATH below to match your system."
    MODULE_PATH="/path/to/GDCGenomicsQC/envs"
    module use "$MODULE_PATH"
    module load gdcgenomicsMSI
    conda activate snakemake
    ;;
  local)
    echo "Setting up local environment..."
    conda activate snakemake
    cd GDCGenomicsQC
    ;;
  *)
    echo "Usage: source setup_env.sh [msi|sandbox|other|local]"
    return 1
    ;;
esac

echo "Environment: $(snakemake --version 2>/dev/null || echo 'snakemake not found')"
