(tutorial_prs)=

# Tutorial: Single-Ancestry Polygenic Risk Score (PRS) Methods

This tutorial covers configuring and running single-ancestry PRS methods using
the GDCGenomicsQC pipeline. Single-ancestry methods train on one ancestry
and apply to the same ancestry.

**Estimated completion time**: 1-2 hours

**Learning objectives**:

1. Configure single-ancestry PRS methods (enable/disable individual methods)
2. Prepare inputs for single-ancestry PRS
3. Run all enabled single-ancestry PRS methods via Snakemake
4. Interpret single-ancestry PRS output files and performance metrics

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
# Replace with your HPC's module path:
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

- Completed [](tutorial_qc_pipeline.md) (QC-filtered genotype data)
- Completed [](tutorial_ancestry_classification.md) (ancestry labels)
- Completed [](tutorial_1kg_assembly.md) (reference data)
- Summary statistics file (from GWAS, format: CHR, BP, SNP, A1, A2, OR, P)
- Phenotype file for target study (format: IID, pheno)
- Optional: LD reference files for PRScs/PRScsx

---

## PRS Methods Overview

The pipeline supports 5 single-ancestry PRS methods (train on one ancestry, apply to same ancestry). Currently, only **PRSice2** and **LDPred2** are fully functional. The roadmap includes enabling the remaining methods in future releases.

```{list-table} Single-Ancestry PRS Methods Status
:widths: 25 25 50
:header-rows: 1

* - Config Key
  - Method Name
  - Status
* - ``single_ct``
  - CT-SLeB
  - Roadmap (not yet functional)
* - ``single_prsice``
  - PRSice2
  - Working
* - ``single_prscs``
  - PRScs
  - Roadmap (not yet functional)
* - ``single_ldpred2``
  - LDPred2
  - Working
* - ``single_lassosum2``
  - lassosum2
  - Roadmap (not yet functional; disabled in config due to missing ``caret`` package in container)
```

All methods are disabled by default. Enable the ones you want to run in the config (only ``single_prsice`` and ``single_ldpred2`` are recommended currently).

For multi-ancestry PRS methods, see [](tutorial_prs_multi.md).

---

### Required Input Files

```{list-table} PRS Input Files
:widths: 35 65
:header-rows: 1

* - Input File
  - Description
* - ``OUT_DIR/{ANC}/f1.b38.f2.pgen``
  - QC-filtered, ancestry-subsetted genotypes
* - ``summary_statistics.tsv``
  - GWAS summary stats (user-provided)
* - ``target_phenotype.tsv``
  - Phenotype for PRS validation (IID, pheno)
* - ``ancestry_labels.tsv``
  - Ancestry labels (from [](tutorial_ancestry_classification.md))
```

**Config Parameters for Single-Ancestry PRS:**

```yaml
# PRS output goes to OUT_DIR/prs/ by default
prsMethods:
  resource_dir: "/path/to/prs_resources"  # Optional, defaults to ../prs_resources
  # Single ancestry methods
  single_ct:
    enabled: false
  single_prsice:
    enabled: true
  single_prscs:
    enabled: false
  single_ldpred2:
    enabled: true
  single_lassosum2:
    enabled: false

conda-frontend: mamba
```

**See also:** [](tutorial_qc_pipeline.md) for genotype prep, [](tutorial_ancestry_classification.md) for ancestry labels.

---

## Lab Exercise: Running All Enabled PRS Methods

### Step 1: Create Configuration File

Create a config that enables multiple PRS methods:

```bash
mkdir -p ~/prs_lab
cd ~/prs_lab
cat > config_prs.yaml << 'EOF'
    INPUT: "/path/to/data/chr{CHR}.vcf.gz"
    REF: "/path/to/reference/storage"
    OUT_DIR: "/path/to/output/directory"
    local-storage-prefix: "/path/to/.snakemake/storage"

    chromosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]

    ancestry:
        ancestry_file: "/path/to/ancestry_labels.tsv"  # Or use model: "pca" for classification

    prsMethods:
        resource_dir: "/path/to/prs_resources"
        # Enable single-ancestry methods
        single_prsice:
          enabled: true
        single_ldpred2:
          enabled: true
        # Disable unused methods
        single_ct:
          enabled: false
        single_prscs:
          enabled: false
        single_lassosum2:
          enabled: false

    conda-frontend: mamba
EOF
```

Key parameters:
- ``prsMethods.<method>.enabled``: Set to ``true`` to run that method
- ``resource_dir``: Directory for LD references and method resources

### Step 2: Run All Enabled Single-Ancestry PRS Methods

The ``runAllEnabledPRS`` target runs all single-ancestry methods marked ``enabled: true`` in config:

::::{tab-set}
:::{tab-item} MSI HPC

```bash
cd GDCGenomicsQC/workflow
gdcgenomicsqc --configfile ../config_prs.yaml runAllEnabledPRS -j 4
```

:::
:::{tab-item} Sandbox

```bash
cd GDCGenomicsQC/workflow
gdcgenomicsqc --configfile ../config_prs.yaml runAllEnabledPRS -j 4
```

:::
:::{tab-item} Other HPCs

```bash
cd GDCGenomicsQC/workflow
gdcgenomicsqc --configfile ../config_prs.yaml runAllEnabledPRS -j 4
```

:::
:::{tab-item} Local Snakemake

```bash
cd GDCGenomicsQC/workflow
snakemake --profile=../profiles/hpc \
    --configfile ../config_prs.yaml \
    runAllEnabledPRS \
    -j 4
```

:::
::::

This will:
1. Prepare PRS resources (LD directories, reference links)
2. Run all enabled single-ancestry methods in parallel
3. Create ``prs_all_completed.done`` when all methods finish

### Step 3: Run Individual Methods (Optional)

To run a single method instead of all enabled:

```bash
# Run only PRSice2
gdcgenomicsqc --configfile ../config_prs.yaml runSingleAncestryPRSice -j 4

# Run only CT-SLeB
gdcgenomicsqc --configfile ../config_prs.yaml runSingleAncestryCT -j 4
```

---

## Interpreting Pipeline Outputs

### Output Directory Structure

PRS outputs are in ``PRS_OUT_DIR/method_runs/``:

```text
prs_output/
├── method_runs/
│   ├── single_ct/
│   │   ├── prs_scores.tsv
│   │   └── performance_metrics.txt
│   ├── single_prsice/
│   │   ├── prsice_summary.csv
│   │   └── best_pRS.prs
│   ├── single_prscs/
│   │   └── prs_scores.tsv
│   ├── single_ldpred2/
│   │   └── prs_scores.tsv
│   └── single_lassosum2/
│       └── prs_scores.tsv
└── prs_all_completed.done
```

### Performance Metrics

**File**: ``method_runs/{method}/performance_metrics.txt``

Sample output:

```{list-table}
:header-rows: 1

* - Method
  - R²
  - AUC
  - p-value
* - single_ct
  - 0.12
  - 0.68
  - 2.3e-5
* - single_prsice
  - 0.15
  - 0.71
  - 1.1e-6
* - single_prscs
  - 0.14
  - 0.69
  - 3.2e-6
```

**Key metrics**:
- ``R²``: Variance in phenotype explained by PRS
- ``AUC``: Area Under ROC Curve (for binary traits)
- ``p-value``: Significance of PRS-phenotype association

---

## Exploration Exercises

Vary these parameters to compare methods:

1. **Method comparison**: Enable all single-ancestry methods and compare R²/AUC across CT-SLeB, PRSice2, PRScs, LDPred2, and lassosum2

2. **P-value threshold**: For methods like PRSice2, test different GWAS p-value cutoffs

3. **Ancestry subsetting**: Run PRS only on EUR vs AFR samples (modify ancestry config)

4. **Resource tuning**: Adjust threads/memory for memory-intensive methods (LDPred2, lassosum2)

---

## Discussion Points

1. **Method performance**: Which single-ancestry method achieves the highest R²/AUC for your trait? Are results consistent with published benchmarks?

2. **LD reference bias**: How do PRS results change with different LD reference panels? What are the implications for underrepresented ancestries?

3. **Computational tradeoffs**: Which methods are fastest? Which require the most memory? How does this affect HPC resource allocation?

4. **P-value thresholding**: How does the optimal p-value threshold vary across methods? What does this reveal about their underlying algorithms?

For theoretical foundations of single-ancestry PRS methods—including LD clumping, p-value thresholding, and Bayesian approaches—refer to accompanying lecture materials.

---

## Next Steps

After completing this tutorial, you have:

- Configured and run single-ancestry PRS methods
- Compared performance across CT-SLeB, PRSice2, PRScs, LDPred2, and lassosum2
- Interpreted PRS output metrics

**Further analyses to consider:**

- Validate PRS in independent holdout samples
- Compare PRS performance across ancestry groups
- For multi-ancestry methods, see [](tutorial_prs_multi.md)

**See also:**

- [](installation.md) - Software setup
- [](tutorial_qc_pipeline.md) - Genotype QC preprocessing
- [](tutorial_ancestry_classification.md) - Ancestry labels for PRS stratification
- [](tutorial_prs_multi.md) - Multi-ancestry PRS methods
- [](genomics.md) - Technical details on PRS methodology

## External Data Inputs

Instead of using simulated phenotypes from the pipeline, you can point to your own
genotype data, GWAS summary statistics, and phenotypes by configuring ``prsPipeline.external``:

```yaml
prsPipeline:
    external:
        target_bed: "/path/to/genotypes"
        target_sumstats: "/path/to/gwas.tsv"
        target_pheno: "/path/to/pheno.tsv"
```

If running multi-ancestry methods, also specify the second ancestry:

```yaml
prsPipeline:
    external:
        target_bed: "/path/to/anc1_genotypes"
        target_sumstats: "/path/to/anc1_gwas.tsv"
        target_pheno: "/path/to/anc1_pheno.tsv"
        anc2_bed: "/path/to/anc2_genotypes"
        training_sumstats: "/path/to/anc2_gwas.tsv"
        training_pheno: "/path/to/anc2_pheno.tsv"
```

**Pipeline-generated files:** If you've run the QC pipeline first, typical paths are:

- Genotypes: ``OUT_DIR/{ANC}/f1.b38.f2.pgen`` (after Standard QC)
- For GRM: ``OUT_DIR/{ANC}/unrelated.grm.bin`` (from Relatedness step)

**GWAS format:** Tab-separated with columns: SNP, A1, A2, BETA, SE, P

**Phenotype format:** Tab-separated with columns: FID, IID, <phenotype>
