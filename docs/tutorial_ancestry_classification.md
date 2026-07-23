(tutorial_ancestry)=

# Tutorial: Ancestry Classification in Practice

This tutorial provides hands-on experience running the ancestry classification
pipeline in GDCGenomicsQC. For the theoretical background on dimension
reduction methods and classification techniques, see the accompanying lecture
slides.

**Estimated completion time**: 30-45 minutes

**Learning objectives**:

1. Run the ancestry classification pipeline using Snakemake
2. Configure different models and thresholds
3. Interpret pipeline outputs
4. Apply ancestry-specific subsetting

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
conda activate snakemake
```

Verify installation:

```bash
snakemake --version
```

```{note}
**You do NOT need to clone the repository.** The pipeline is pre-installed
via the ``gdcgenomicsMSI`` module. Just create your config file and run.
```

:::
:::{tab-item} Sandbox

If you're using the Sandbox environment:

```bash
module use /scratch.global/GDC/GDCGenomicsQC/envs
module load gdcgenomicsSandbox
conda activate snakemake
```

Verify installation:

```bash
snakemake --version
```

```{note}
**You do NOT need to clone the repository.** The pipeline is pre-installed
via the ``gdcgenomicsSandbox`` module. Just create your config file and run.
```

:::
:::{tab-item} Other HPCs

If your HPC has the GDC module pre-configured:

```bash
# Replace with your HPC's module path:
module use /path/to/GDCGenomicsQC/envs
module load gdcgenomicsMSI
conda activate snakemake
```

Verify installation:

```bash
cd GDCGenomicsQC
snakemake --version
```

:::
:::{tab-item} Local Snakemake

If you're using your own Snakemake installation:

```bash
conda activate snakemake
cd GDCGenomicsQC
```

Verify installation:

```bash
snakemake --version
```

:::
::::

**Data Requirements:**

- Reference panel with population labels (see [](tutorial_1kg_assembly.md))
- QC-filtered genotype data (see [](tutorial_qc_pipeline.md))

(dag-visualization)=

### DAG Visualization

The pipeline DAG up to the ``run_pca`` rule shows the workflow for preparing
the PCA reference and running ancestry classification:

```{mermaid} dag_pca.mmd

```

The file graph shows how input and output files flow through each rule:

```{mermaid} filegraph_pca.mmd

```

---

### Required Input Files

This step requires the following input files:

```{list-table} Ancestry Classification Input Files
:widths: 35 65
:header-rows: 1

* - Input File
  - Description
* - ``INPUT: "chr{CHR}.vcf.gz"`` (or .bed/.pgen)
  - Per-chromosome genotype data (QC-filtered recommended)
* - ``REF/1000G_highcoverage/population.txt``
  - Reference panel with population labels (IID, pop, superpop columns)
* - ``REF/1000G_highcoverage/1000G_highCoveragephased.pruned.pgen``
  - LD-pruned, unrelated reference genotypes for PCA projection
* - ``OUT_DIR/full/f1.pgen`` (or ``f1_{CHR}.pgen``)
  - Initial QC-filtered sample genotypes
```

**Input from Previous Steps:**

The ancestry classification pipeline depends on:

1. **QC Pipeline** (tutorial_qc_pipeline): Produces filtered genotype files
2. **Reference Assembly** (tutorial_1kg_assembly): Provides reference panel

**Config Parameters for Ancestry:**

```yaml
ancestry:
    threshold: 0.8  # Minimum posterior probability for classification
    model: "pca"    # Options: pca, umap, rfmix (vae not yet implemented)
    pca_estimation: "projection"  # "projection" or "joint" — how PCA is computed
    # Optional: reported_race: "/path/to/reported_race.tsv"

# Optional: subset samples/variants at the very start of the pipeline
keep_samples: "/path/to/sample_iids.txt"       # One IID per line
keep_variants: "/path/to/variant_ids.txt"      # One variant ID per line
remove_samples: "/path/to/remove_iids.txt"     # One IID per line — removed before keep files
exclude_variants: "/path/to/exclude_vars.txt"  # One variant ID per line — excluded before keep_variants

# QC missingness thresholds (optional, default values shown)
initial_variant_missingness: 0.1     # Initial --geno threshold (convertPlink per-chromosome step)
final_variant_missingness: 0.02      # Final --geno threshold (initialFilter.sh)
initial_subject_missingness: 0.1     # Initial --mind threshold (initialFilter.sh)
final_subject_missingness: 0.02      # Final --mind threshold (initialFilter.sh)

# HWE sample-size scaling factor (optional, default null = fixed threshold)
# hwe_k: 0.001   # Greer et al. (2024) recommends 0.001 for large studies

INPUT: "/path/to/data/chr{CHR}.vcf.gz"
OUT_DIR: "/path/to/output"
REF: "/path/to/reference"
local-storage-prefix: "/path/to/.snakemake/storage"

chromosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
```

**See also:** [](tutorial_qc_pipeline.md) for QC preprocessing, [](tutorial_1kg_assembly.md) for reference data.

## Lab Exercise: Running Ancestry Classification

### Step 1: Create Configuration File

For this tutorial, we will classify ancestry using an admix-free resampling
of the 1000 Genomes reference panel. This reference was generated using the
CT-Sleb tool and is available on Harvard Dataverse. The resampling ensures
that only genetically unrelated, ancestry-appropriate samples are included,
providing a cleaner reference for classification.

Create a configuration file for ancestry classification:

```bash
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
    model: "pca"  # Options: pca, umap, rfmix (vae not yet implemented)
    pca_estimation: "projection"  # "projection" or "joint"

# Optional: subset samples/variants before ancestry classification
# keep_samples: "/path/to/sample_iids.txt"
# keep_variants: "/path/to/variant_ids.txt"
# remove_samples: "/path/to/remove_iids.txt"
# exclude_variants: "/path/to/exclude_vars.txt"

# Optional: QC missingness thresholds (defaults shown)
# initial_variant_missingness: 0.1
# final_variant_missingness: 0.02
# initial_subject_missingness: 0.1
# final_subject_missingness: 0.02
# hwe_k: null    # set to 0.001 for Greer et al. (2024) recommendation

relatedness:
    method: "king"  # "0" for none, "king" or "primus" for removal
    king_cutoff: 0.0884

# Internal PCA (optional)
internalPCA:
    method: "plink2"  # "plink2", "pcair", or "both"
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
```

Key parameters:

- ``threshold``: Minimum posterior probability for confident classification (default: 0.8)
- ``model``: Embedding used for classification—``pca``, ``umap``, or ``rfmix``
  (Note: VAE is not yet implemented)
- ``pca_estimation``: How PCA components are computed:
  - ``"projection"`` (default): PCA on the 1000G reference panel only, then projects study samples onto those PCs. Fast, reference-consistent.
  - ``"joint"``: Merges study and reference genotypes, computes PCA jointly, then splits by population. Better for capturing study-specific variation but slower.
- ``keep_samples``: Path to a file with sample IIDs (one per line) to subset data at the start of the pipeline. Applied on top of ancestry-specific keep files.
- ``keep_variants``: Path to a file with variant IDs (one per line) to subset variants at the start of the pipeline.
- ``remove_samples``: Path to a file with sample IIDs (one per line) to remove from the data at the start. These samples are excluded before ancestry-specific keep files are applied.
- ``exclude_variants``: Path to a file with variant IDs (one per line) to exclude from the data at the start. Applied before ``keep_variants`` (exclusion first, then extraction). IDs must match the variant ID format in the input data (e.g., ``chr:pos:ref:alt``).
- ``initial_variant_missingness``: Initial --geno threshold (default: 0.1). Applied per-chromosome in ``convertPlink`` before allele alignment. Removes variants with >10% missing genotypes.
- ``final_variant_missingness``: Final --geno threshold (default: 0.02). Applied in ``initialFilter.sh`` after allele alignment. More stringent, removes variants with >2% missing genotypes.
- ``initial_subject_missingness``: Initial --mind threshold (default: 0.1). Applied in ``initialFilter.sh`` to remove samples with >10% missing genotypes.
- ``final_subject_missingness``: Final --mind threshold (default: 0.02). Applied in ``initialFilter.sh`` after variant missingness filter. Removes samples with >2% missing genotypes.
- ``hwe_k``: Sample-size scaling factor k for --hwe in ``applyStandardQualityControl`` (default: null, meaning k=0/fixed threshold). The effective HWE p-value becomes p × 10^(−n×k). Greer et al. (2024) recommends k=0.001 for large studies to avoid discarding genuine associations from overpowered HWE tests.

### Step 2: Run Classification Pipeline

Before classification, the pipeline runs KING to identify related samples.
For this tutorial using simulated data, you should not find any related
individuals (KING kinship coefficient ≈ 0), which serves as a good sanity
check that the simulated data is properly independent.

::::{tab-set}
:::{tab-item} MSI HPC

```bash
cd GDCGenomicsQC/workflow
gdcgenomicsqc --configfile ../config_ancestry.yaml classifyAncestry -j 10
```

:::
:::{tab-item} Sandbox

```bash
cd GDCGenomicsQC/workflow
gdcgenomicsqc --configfile ../config_ancestry.yaml classifyAncestry -j 10
```

:::
:::{tab-item} Other HPCs

```bash
cd GDCGenomicsQC/workflow
gdcgenomicsqc --configfile ../config_ancestry.yaml classifyAncestry -j 10
```

:::
:::{tab-item} Local Snakemake

```bash
cd GDCGenomicsQC/workflow
snakemake --profile=../profiles/hpc \
    --configfile ../config_ancestry.yaml \
    classifyAncestry \
    -j 10
```

:::
::::

This trains Random Forest models on reference coordinates and predicts ancestry
probabilities for your samples.

### Step 3: Compare Models (PCA vs UMAP)

Modify ``model`` in your config to compare embeddings:

- **PCA** (default): Linear projection, strongest baseline
- **UMAP**: Nonlinear, good for visualization
- **VAE**: Not yet implemented

::::{tab-set}
:::{tab-item} MSI HPC

First edit your config to set ``model: "umap"``, then:

```bash
gdcgenomicsqc --configfile ../config_ancestry.yaml classifyAncestry
```

:::
:::{tab-item} Sandbox

First edit your config to set ``model: "umap"``, then:

```bash
gdcgenomicsqc --configfile ../config_ancestry.yaml classifyAncestry
```

:::
:::{tab-item} Other HPCs

First edit your config to set ``model: "umap"``, then:

```bash
gdcgenomicsqc --configfile ../config_ancestry.yaml classifyAncestry
```

:::
:::{tab-item} Local Snakemake

First edit your config to set ``model: "umap"``, then:

```bash
snakemake --profile=../profiles/hpc \
    --configfile ../config_ancestry.yaml \
    classifyAncestry
```

:::
::::

### Step 4: Ancestry-Specific Subsetting

The pipeline creates keep files for each predicted ancestry:

::::{tab-set}
:::{tab-item} MSI HPC

```bash
gdcgenomicsqc --configfile ../config_ancestry.yaml convertNfilt/CHR=20/subset=EUR
```

:::
:::{tab-item} Sandbox

```bash
gdcgenomicsqc --configfile ../config_ancestry.yaml convertNfilt/CHR=20/subset=EUR
```

:::
:::{tab-item} Other HPCs

```bash
gdcgenomicsqc --configfile ../config_ancestry.yaml convertNfilt/CHR=20/subset=EUR
```

:::
:::{tab-item} Local Snakemake

```bash
snakemake --profile=../profiles/hpc \
    --configfile ../config_ancestry.yaml \
    convertNfilt/CHR=20/subset=EUR
```

:::
::::

Available subsets are dynamically determined from classification results.

---

## Interpreting Pipeline Outputs

### Classification Probabilities

**File**: ``01-globalAncestry/classificationProbabilities.tsv``

Sample output:

```{list-table}
:header-rows: 1

* - IID
  - pca_AFR
  - pca_AMR
  - pca_EUR
  - pca_SAS
* - Sample1
  - 0.95
  - 0.02
  - 0.02
  - 0.01
* - Sample2
  - 0.05
  - 0.10
  - 0.83
  - 0.02
* - Sample3
  - 0.40
  - 0.30
  - 0.15
  - 0.15
```

### Classifications

**File**: ``01-globalAncestry/ancestry_classifications.tsv``

```{list-table}
:header-rows: 1

* - IID
  - pca_predicted
  - pca_confidence
* - Sample1
  - AFR
  - 0.95
* - Sample2
  - EUR
  - 0.83
* - Sample3
  - uncertain
  - 0.40
```

Samples below threshold are labeled "uncertain" or grouped as "Other".

### Keep Files

PLINK-style files for ancestry-specific analyses:

- ``keep_AFR.txt``, ``keep_EUR.txt``, etc.
- ``keep_Other.txt`` (below threshold)

Each keep file has two columns (``FID IID``) as required by ``plink2 --keep``.
If the input data has no FID column, ``FID`` is set equal to ``IID``.

### Visualizations

**Stacked Area Plot**: ``classificationProbability_stacked_pca.svg``

- X-axis: Samples sorted by ancestry proportions
- Y-axis: Stacked classification probabilities
- Identifies homogeneous and admixed individuals

**Classification Space**: ``images/ancestry_classification_space.svg``

- A visualization of classification using PCA
- Samples in PC space with reference density contours
- Color indicates predicted ancestry

### Creating a Confusion Matrix with Reported Race Labels

To evaluate classification performance, you can compare predicted ancestry
labels against reported race/ethnicity data. Provide a tab-separated file
with sample IDs and reported labels:

**Input format** (``reported_race.tsv``):

```{list-table}
:header-rows: 1

* - IID
  - reported
* - Sample1
  - AFR
* - Sample2
  - EUR
* - Sample3
  - unknown
```

To generate the confusion matrix, add to your config:

```bash
ancestry:
    threshold: 0.8
    model: "pca"
    reported_race: "/path/to/reported_race.tsv"
```

The pipeline will output:

- ``ancestry_confusion_matrix.tsv``: Contingency table of predicted vs. reported
- ``ancestry_confusion_matrix.svg``: Heatmap visualization

**Interpretation notes**:

- Self-reported race is a social construct, not a genetic one—expect imperfect
  concordance due to genetic ancestry not aligning with social categorization
- Admixed individuals may not map cleanly to discrete categories
- Discrepancies can reveal both classification errors and limitations of
  self-reported labels

### Using a Provided Ancestry Classification File

If you already have ancestry labels for your samples (e.g., from a previous
analysis, clinical database, or external classifier), you can bypass the
pipeline's ancestry prediction entirely by providing a tab-separated file.

**When to use this:**

- You have existing ancestry labels you trust
- You want faster pipeline execution (skips PCA/UMAP/RFMix)
- You need specific ancestry labels not supported by the default classifier

**Input format** (``ancestry_labels.tsv``):

```{list-table}
:header-rows: 1

* - IID
  - ancestry
* - Sample1
  - AFR
* - Sample2
  - EUR
* - Sample3
  - EUR
```

The file should be:

- Tab-separated (or specify a different separator with ``ancestry_file_sep``)
- Two columns: IID (sample ID), ancestry label
- No header row (or specify a column name with ``ancestry_file_col``)
- One line per sample

To use your labels, add to your config:

```bash
ancestry:
    threshold: 0.8
    model: "pca"
    ancestry_file: "/path/to/ancestry_labels.tsv"
    # Optional: ancestry_file_sep: " "         # for space-separated files
    # Optional: ancestry_file_col: "Superpopulation"  # column name in header
```

**How the bypass works:**

1. The pipeline reads your file and extracts unique ancestry labels
2. Creates ``keep_{ancestry}.txt`` files in ``01-globalAncestry/`` (same as predicted)
3. Skips ancestry prediction rules (PCA, UMAP, RFMix outputs are not required)
4. Branches ancestry-specific QC using your provided labels

**Behavior:**

- Samples NOT in your file are excluded from ancestry-specific QC
- They remain in the "full" dataset for non-stratified analyses
- The ``phenotypeSimulation.ancestries`` config must match labels in your file

**Example complete config:**

```bash
INPUT: "/path/to/data/chr{CHR}.vcf.gz"
OUT_DIR: "/path/to/output/directory"
REF: "/path/to/reference/data"
local-storage-prefix: "/path/to/.snakemake/storage"

chromosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]

ancestry:
    model: "pca"
    ancestry_file: "/path/to/ancestry_labels.tsv"

phenotypeSimulation:
    ancestries: ["AFR", "EUR"]  # Must match labels in your file
```

This enables rapid iteration when you already have ancestry assignments.

---

## Discussion Points

These questions extend the practical exercise into deeper methodological considerations:

1. **Model comparison**: How do posterior probability distributions differ between
   PCA and UMAP? Does this align with the simulation findings that PCA
   remains the strongest baseline? (VAE not yet available for comparison)

2. **Threshold selection**: What happens to the number of "uncertain" classifications
   as you vary the threshold from 0.6 to 0.95? How does this affect downstream
   sample sizes?

3. **Admixed samples**: Examine samples with mixed ancestry proportions in the
   stacked area plot. Should these be forced into discrete categories, or would
   soft probabilities be more appropriate for covariate adjustment?

4. **Reference panel bias**: How do classifications change if your target
   population differs from the reference panel? What are the implications for
   fairness and validity?

5. **Classification vs. covariates**: For GWAS adjustment, compare results using
   hard ancestry labels versus PCs as continuous covariates. Which approach is
   more appropriate and why?

6. **Confusion and error**: Which ancestry pairs are most frequently confused
   in your data? Is this consistent with the simulation results showing PCA as
   nearly perfect on pure-like samples?

7. **Uncertainty quantification**: The pipeline provides probability estimates.
   How should these be incorporated into downstream analyses? Should low-confidence
   samples be excluded or modeled differently?

For the theoretical foundations behind these methods—including PCA decomposition,
Random Forest ensemble learning, and evaluation metrics—refer to the accompanying
lecture materials.

---

## Next Steps

After completing this tutorial, you can:

- [](tutorial_heritability.md) - Estimate heritability using ancestry-classified samples
- Return to [](tutorial_qc_pipeline.md) - Run ancestry-specific QC using the keep files

## **Lab Materials**

- [Global Ancestry Visualization Lab (Quarto)](labs/lab03_global_ancestry_visualization.qmd) - Interactive R notebook for visualizing global ancestry outputs

**The ancestry classification outputs enable:**

- Ancestry-specific QC filters (``EUR/f1.b38.f2.pgen``, ``AFR/f1.b38.f2.pgen``, etc.)
- Per-ancestry heritability estimation
- Stratified GWAS analyses

**See also:**

- [](installation.md) - Software setup (if not already done)
- [](usage.md) - Running the full pipeline
- [](genomics.md) - Technical details on ancestry methods
