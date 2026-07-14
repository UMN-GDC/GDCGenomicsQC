# Welcome to the GDC Documentation!

This project provides tools and pipelines for interacting with the
Genomic Data Commons (GDC). It is designed to run efficiently on
high-performance clusters like the UMN MSI Agate cluster.

````{toctree}
:maxdepth: 2
:caption: Contents:

installation
new_hpc_setup
genomics
usage
tutorial_qc_pipeline
tutorial_ancestry_classification
tutorial_1kg_assembly
tutorial_heritability
tutorial_phenotype_simulation
dev_environments
````

## Getting Started

**Recommended Learning Pathway:**

1. [](installation.md) - Set up software environment
2. [](usage.md) - Learn how to run the pipeline
3. [](tutorial_1kg_assembly.md) - Download reference data
4. [](tutorial_qc_pipeline.md) - Run quality control
5. [](tutorial_ancestry_classification.md) - Classify ancestry
6. [](tutorial_phenotype_simulation.md) - Simulate phenotypes for method testing
7. [](tutorial_heritability.md) - Estimate heritability (with real phenotypes)

**Quick Setup (MSI/UMN HPC):**

```bash
module use /path/to/GDCGenomicsQC/envs
module load gdcgenomicsMSI
conda activate snakemake
cd GDCGenomicsQC
snakemake --version
```

**See also:** [](usage.md) for detailed instructions on running the pipeline with module load or local snakemake.

## Indices and tables

- {ref}`genindex`
- {ref}`modindex`
- {ref}`search`
