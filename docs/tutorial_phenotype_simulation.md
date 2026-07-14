(tutorial_phenotype_simulation)=

# Tutorial: Phenotype Simulation and Heritability Estimation

This tutorial covers simulating phenotypes across multiple ancestries and
estimating SNP heritability using the GDCGenomicsQC pipeline.

**Estimated completion time**: 2-3 hours

**Learning objectives**:

1. Configure phenotype simulation for multiple ancestry groups
2. Simulate phenotypes with controlled heritability and cross-ancestry genetic correlation
3. Run SNP heritability estimation using PC-relate
4. Compare heritability estimates across ancestries

---

## Prerequisites

**Setup:**

Before starting, ensure you have access to Snakemake and the GDCGenomicsQC workflow.
For detailed installation instructions, see:

- [](installation.md) - Software setup (module, conda, or other methods)
- [](usage.md) - Running the pipeline

::::{tab-set}
:::{tab-item} MSI HPC

If you're using the MSI HPC cluster:

```bash
module use /projects/standard/gdc/public/GDCGenomicsQC/envs
module load gdcgenomicsMSI
```

:::
:::{tab-item} Sandbox

If you're using the Sandbox environment:

```bash
module use /scratch.global/GDC/GDCGenomicsQC/envs
module load gdcgenomicsSandbox
```

:::
:::{tab-item} Other HPCs

If your HPC has the GDC module pre-configured:

```bash
module use /path/to/GDCGenomicsQC/envs
module load gdcgenomicsMSI
```

:::
:::{tab-item} Local Snakemake

If you're using your own Snakemake installation:

```bash
conda activate snakemake
cd GDCGenomicsQC
```

:::
::::

**Data Requirements:**

- Completed [](tutorial_1kg_assembly.md) (reference data assembled)
- QC-filtered genotype data (see [](tutorial_qc_pipeline.md))
- Ancestry labels (either from [](tutorial_ancestry_classification.md) or provided file)

---

## Phenotype Simulation Overview

The phenotype simulation feature generates synthetic phenotypes with
controlled genetic architecture for testing heritability estimation methods.

**Key capabilities:**

- Simulate phenotypes with specified SNP heritability (h²)
- Control cross-ancestry genetic correlation via an N×N correlation matrix
- Generate multiple independent simulations per configuration
- Support N ancestries (not limited to 2-way bivariate models)
- Multiple named simulations per pipeline run

---

### Required Input Files

This step requires the following input files:

```{list-table} Phenotype Simulation Input Files
:widths: 35 65
:header-rows: 1

* - Input File
  - Description
* - ``OUT_DIR/{ANC}/f1.pgen``
  - Sample genotypes per ancestry (from QC pipeline)
* - ``ancestry_file`` (optional)
  - User-provided ancestry labels (see [](tutorial_ancestry_classification.md))
```

**Two ways to provide ancestry:**

1. **Run ancestry classification** - Use [](tutorial_ancestry_classification.md) to predict labels
2. **Provide your own** - Use ``ancestry_file`` config option (faster)

**Simulation Parameters:**

```yaml
phenotypeSimulation:
    enabled: true              # Set to true to run simulation
    ancestries: ["AFR", "EUR"] # Ancestry groups to simulate
    simulations:
        - name: "sim1"
          # N×N genetic correlation matrix (size = len(ancestries))
          corr_matrix: [[1.0, 0.8], [0.8, 1.0]]
          n_sims: 10           # Number of phenotype simulations
          heritability: 0.4    # SNP heritability (h²)
          maf: 0.05            # Minor allele frequency threshold
          seed: 42             # Random seed for reproducibility
          skip_thinning: true  # Skip SNP thinning
```

**Heritability Estimation on Simulated Data:**

```yaml
snpHerit:
    # No phenotype needed — estimateSnpHeritabilitySimulated
    # uses the simulated phenotype automatically
    method: "AdjHE"   # Estimation method
    npc: 10           # Number of PCs to include
```

**Output Files:**

```{list-table} Simulation Output Files
:widths: 40 60
:header-rows: 1

* - File
  - Description
* - ``{ancestry}/simulations/{sim_name}/simulated.bed``
  - Simulated genotype PLINK files (per ancestry)
* - ``{ancestry}/simulations/{sim_name}/simulated_pheno1.pheno``
  - Simulated phenotype file (per ancestry)
* - ``{ancestry}/simulations/{sim_name}/herit.csv``
  - Heritability estimates from MASH (per ancestry)
```

---

## Lab Exercise: Phenotype Simulation

### Option A: Using Provided Ancestry Labels (Recommended for Testing)

This approach skips ancestry prediction and is faster for testing.

#### Step 1: Create ancestry labels file

Create a tab-separated file with sample IDs and ancestry labels:

```bash
mkdir -p ~/sim_lab
cd ~/sim_lab
cat > ancestry_labels.tsv << 'EOF'
sample1	AFR
sample2	AFR
sample3	EUR
sample4	EUR
EOF
```

#### Step 2: Create configuration file

```bash
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
    simulations:
        - name: "sim1"
          corr_matrix: [[1.0, 0.8], [0.8, 1.0]]
          heritability: 0.4
          n_sims: 10
          skip_thinning: true

snpHerit:
    method: "AdjHE"
    npc: 10

conda-frontend: mamba
EOF
```

#### Step 3: Run simulation

::::{tab-set}
:::{tab-item} MSI HPC

```bash
cd GDCGenomicsQC/workflow
gdcgenomicsqc --configfile ../config_simulation.yaml run_simulatePhenotype -j 4
```

:::
:::{tab-item} Sandbox

```bash
cd GDCGenomicsQC/workflow
gdcgenomicsqc --configfile ../config_simulation.yaml run_simulatePhenotype -j 4
```

:::
:::{tab-item} Local Snakemake

```bash
cd GDCGenomicsQC/workflow
snakemake --profile=../profiles/hpc \
    --configfile ../config_simulation.yaml \
    run_simulatePhenotype \
    -j 4
```

:::
::::

Output directory structure:

```text
AFR/simulations/sim1/
├── simulated.bed
├── simulated.bim
├── simulated.fam
├── simulated_pheno1.pheno
├── simulated.grm.bin
├── simulated.grm.id
├── simulated.grm.N.bin
├── simulated.eigenvec
└── herit.csv
EUR/simulations/sim1/
├── simulated.bed
├── simulated.bim
├── simulated.fam
├── simulated_pheno1.pheno
├── simulated.grm.bin
├── simulated.grm.id
├── simulated.grm.N.bin
├── simulated.eigenvec
└── herit.csv
```

### Option B: Using Predicted Ancestry Labels

If you want to use the full ancestry classification pipeline:

#### Step 1: First run ancestry classification

See [](tutorial_ancestry_classification.md) to run the full classification.

#### Step 2: Update config to use predicted ancestries

Remove the ``ancestry_file`` line and the pipeline will use predicted labels:

```bash
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
    simulations:
        - name: "sim1"
          corr_matrix: [[1.0, 0.8], [0.8, 1.0]]
          heritability: 0.4
          n_sims: 10

snpHerit:
    method: "AdjHE"
    npc: 10

conda-frontend: mamba
EOF
```

#### Step 3: Run simulation

```bash
cd GDCGenomicsQC/workflow
gdcgenomicsqc --configfile ../config_simulation_predicted.yaml run_simulatePhenotype -j 4
```

This will run the full ancestry classification pipeline first, then proceed
to phenotype simulation.

## Simulation Parameters

```{list-table}
:widths: 25 25 50
:header-rows: 1

* - Parameter
  - Default
  - Description
* - ``ancestries``
  - [AFR, EUR]
  - Ancestry groups to simulate
* - ``simulations[].name``
  - required
  - Simulation name (subdirectory name)
* - ``simulations[].corr_matrix``
  - required
  - N×N genetic correlation matrix
* - ``simulations[].n_sims``
  - 10
  - Number of phenotype simulations
* - ``simulations[].heritability``
  - 0.4
  - SNP heritability (h²) for each ancestry
* - ``simulations[].maf``
  - 0.05
  - Minor allele frequency threshold
* - ``simulations[].seed``
  - 42
  - Random seed for reproducibility
* - ``simulations[].skip_thinning``
  - true
  - Skip SNP thinning (faster)
* - ``simulations[].thin_count_snps``
  - 1000000
  - SNPs to thin to (if not skipping)
* - ``simulations[].thin_count_inds``
  - 10000
  - Individuals to thin to (if not skipping)
```

---

## Interpreting Results

### Simulation Results

**Per-ancestry file**: ``{ancestry}/simulations/{sim_name}/herit.csv``

Sample output:

```{list-table}
:header-rows: 1

* - Ancestry
  - h2
  - SE
* - AFR
  - 0.38
  - 0.05
* - EUR
  - 0.42
  - 0.04
```

**To run heritability estimation on simulated data:**

```bash
gdcgenomicsqc --configfile config_simulation.yaml run_snpHeritSimulated -j 4
```

This invokes the ``estimateSnpHeritabilitySimulated`` rule for each
ancestry × simulation combination.

### Expected Results

Given simulation parameters:

- True h² = 0.4 (specified)
- Expected estimates: 0.35-0.45 (within sampling error)
- Cross-ancestry ρ = 0.8 (from corr_matrix)

Differences between ancestries may reflect:

- LD score differences
- Sample size variation
- Genetic architecture heterogeneity

### Comparison Across Simulations

With multiple simulations (``n_sims: 10``), you can analyze:

- Distribution of heritability estimates
- Standard error of estimates
- Bias in estimation method

---

## Exploration Exercises

Vary these parameters to understand the methods:

1. **Heritability**: Test h² = 0.1, 0.3, 0.5, 0.7
   - How does estimation accuracy change?

2. **Cross-ancestry correlation**: Test different values in corr_matrix
   - What happens when ρ = 1 (identical genetic architecture)?

3. **Sample size**: Vary ``thin_count_inds``
   - How does precision improve with more samples?

4. **SNP density**: Vary ``thin_count_snps``
   - Effect of SNP count on heritability estimates

5. **MAF threshold**: Test maf = 0.01, 0.05, 0.10
   - Impact of rare variant inclusion

6. **N-ancestry**: Add a third ancestry group
   - Requires a 3×3 corr_matrix

---

## Discussion Points

1. **Estimation bias**: How close are the estimated h² values to the true
   simulated value (0.4)? What factors contribute to bias?

2. **Method comparison**: Compare AdjHE vs. other methods (GCTA, PredLMM).
   Which is more accurate? More precise?

3. **Cross-ancestry correlation**: When ρ < 1, what does this imply about
   genetic architecture differences? How does this affect meta-analysis?

4. **Sample size effects**: How do standard errors change with sample size?
   Is there a point of diminishing returns?

5. **PC correction**: How many PCs are optimal for controlling population
   structure? What happens with too few or too many?

6. **Heritability heterogeneity**: Why might h² differ between AFR and EUR
   even when simulated with the same true value?

For the theoretical foundations of SNP heritability estimation—including
PC-relate methodology, REML estimation, and genetic correlation—refer to
the accompanying lecture materials.

---

## Next Steps

After completing this tutorial, you have explored:

- Phenotype simulation with controlled heritability
- SNP heritability estimation using PC-relate
- Cross-ancestry genetic correlation
- Using provided ancestry labels for faster pipeline execution

**Further analyses to consider:**

- GWAS on simulated phenotypes with known true effects
- Compare heritability estimates across different ancestry groups
- Test different heritability estimation methods and PC covariates
- Add a third ancestry group with a 3×3 corr_matrix

**See also:**

- [](installation.md) - Software setup
- [](tutorial_ancestry_classification.md) - Using provided ancestry labels
- [](tutorial_qc_pipeline.md) - QC pipeline
- [](genomics.md) - Technical details on heritability methods
