#!/bin/bash
# Phenotype Simulation Pipeline
#
# Simulates phenotypes with controlled heritability and cross-ancestry
# genetic correlation for testing heritability estimation methods.
#
# Usage: bash run_phenotype_sim.sh [msi|sandbox|other|local] [option_a|option_b]
#
# Prerequisites:
# - Run setup_env.sh first (or source it)
# - Ancestry labels available (Option A) or classification completed (Option B)

ENV="${1:-msi}"
OPTION="${2:-option_a}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$SCRIPT_DIR/setup_env.sh" "$ENV"

mkdir -p ~/sim_lab
cd ~/sim_lab

if [ "$OPTION" = "option_a" ]; then
  # Option A: Using provided ancestry labels (faster)
  echo "Option A: Using provided ancestry labels"

  # Create ancestry labels file
  cat > ancestry_labels.tsv << 'EOF'
sample1  AFR
sample2  AFR
sample3  EUR
sample4  EUR
EOF

  # Create config
  cat > config_simulation.yaml << 'EOF'
INPUT: "/path/to/data/chr{CHR}.vcf.gz"
REF: "/path/to/reference/storage"
OUT_DIR: "/path/to/output/directory"
local-storage-prefix: "/path/to/.snakemake/storage"

chromosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]

ancestry:
    ancestry_file: "/path/to/ancestry_labels.tsv"

phenotypeSimulation:
    enabled: true
    ancestries: ["AFR", "EUR"]
    simulations_dir: "/path/to/simulations"
    n_sims: 10
    heritability: 0.4
    rho: 0.8
    maf: 0.05
    seed: 42
    skip_thinning: true

snpHerit:
    method: "AdjHE"
    npc: 10

conda-frontend: mamba
EOF

  echo "Config written to config_simulation.yaml"

  cd GDCGenomicsQC/workflow
  echo "Running simulation..."
  gdcgenomicsqc --configfile ../config_simulation.yaml simulatePhenotype -j 4

else
  # Option B: Using predicted ancestry from classification
  echo "Option B: Using predicted ancestry from classification"

  cat > config_simulation_predicted.yaml << 'EOF'
INPUT: "/path/to/data/chr{CHR}.vcf.gz"
REF: "/path/to/reference/storage"
OUT_DIR: "/path/to/output/directory"
local-storage-prefix: "/path/to/.snakemake/storage"

chromosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]

ancestry:
    model: "pca"
    threshold: 0.8

phenotypeSimulation:
    ancestries: ["AFR", "EUR"]
    n_sims: 10
    heritability: 0.4
    rho: 0.8
    seed: 42

snpHerit:
    method: "AdjHE"
    npc: 10
    out: "heritability_estimates.txt"

conda-frontend: mamba
EOF

  echo "Config written to config_simulation_predicted.yaml"

  cd GDCGenomicsQC/workflow
  echo "Running simulation..."
  gdcgenomicsqc --configfile ../config_simulation_predicted.yaml simulatePhenotype -j 4
fi

echo "Phenotype simulation complete."
echo "Outputs in: simulations/{ANC1}_{ANC2}/"
