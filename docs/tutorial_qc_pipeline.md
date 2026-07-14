(tutorial_qc)=

# Tutorial: Quality Control Pipeline in Practice

This tutorial provides hands-on experience running the quality control (QC)
pipeline in GDCGenomicsQC. The pipeline implements a two-stage approach:
Initial QC (sample and variant missingness filtering) followed by Standard QC
(MAF, HWE, heterozygosity, and optional sex checking).

**Estimated completion time**: 20-30 minutes

**Learning objectives**:

1. Run the Initial QC and Standard QC pipelines
2. Interpret output plots and data files
3. Configure QC thresholds for different study designs
4. Apply ancestry-specific QC workflows

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

- Reference data configured (see [](tutorial_1kg_assembly.md))
- Genotype data in VCF, BED, or PGEN format

(dag-visualization-qc)=

### DAG Visualization

The pipeline DAG up to the ``run_initialQC`` rule shows the Initial QC workflow:

```{mermaid} dag_initialQC.mmd

```

The rule graph provides a cleaner view of rule dependencies:

```{mermaid} rulegraph_initialQC.mmd

```

---

### Required Input Files

This step requires the following input files:

```{list-table} QC Pipeline Input Files
:widths: 35 65
:header-rows: 1

* - Input File
  - Description
* - ``INPUT: "chr{CHR}.vcf.gz"`` (or .bed/.pgen)
  - Per-chromosome VCF, BED, or PGEN files with genotype data
* - ``REF/1000G_highcoverage/population.txt``
  - Reference panel population labels (for ancestry QC subsets)
* - ``REF/Homo_sapiens.GRCh38.dna.primary_assembly.fa``
  - Reference genome FASTA (if using reference allele correction)
```

**Input Formats Supported:**

The pipeline automatically detects format based on file extension:

```{list-table}
:header-rows: 1

* - Format
  - Example Path
* - VCF
  - ``/data/chr{CHR}.vcf.gz``
* - PLINK BED
  - ``/data/chr{CHR}.bed``
* - PLINK PGEN
  - ``/data/chr{CHR}.pgen``
* - Single file
  - ``/data/merged.bed`` (no ``{CHR}``)
```

**Config Parameters for QC:**

```yaml
INPUT: "/path/to/data/chr{CHR}.vcf.gz"  # Per-chromosome VCF
OUT_DIR: "/path/to/output"
REF: "/path/to/reference"
local-storage-prefix: "/path/to/.snakemake/storage"

chromosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]

# QC thresholds
relatedness:
    method: "king"  # "0" for none, "king" for removal
    king_cutoff: 0.0884

SEX_CHECK: true  # Enable/disable sex verification
GRM: true  # Compute genetic relationship matrix
thin: false
```

**See also:** [](usage.md) for configuration options, [](installation.md) for software setup.

---

## Lab Exercise: Running QC Pipeline

### Step 1: Create Configuration File

Create a configuration file for QC:

```bash
mkdir -p ~/qc_lab
cd ~/qc_lab
cat > config_qc.yaml << 'EOF'
INPUT: "/path/to/data/chr{CHR}.vcf.gz"
OUT_DIR: "/path/to/output/directory"
REF: "/path/to/reference/data"
local-storage-prefix: "/path/to/.snakemake/storage"

chromosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]

relatedness:
    method: "king"  # "0" for none, "king" or "primus" for removal
    king_cutoff: 0.0884

SEX_CHECK: true
thin: false
conda-frontend: mamba

# Internal PCA (optional)
internalPCA:
    method: "plink2"  # "plink2", "pcair", or "both"
    npc: 20
EOF
```

Key parameters:

- ``SEX_CHECK``: Enable/disable sex verification (default: true)
- ``relatedness.method``: Relatedness filtering method ("0" for none, "king" or "primus" for removal)
- ``internalPCA.method``: PCA method ("plink2", "pcair", or "both")

### Step 2: Run Initial QC

::::{tab-set}
:::{tab-item} MSI HPC

```bash
cd GDCGenomicsQC/workflow
gdcgenomicsqc --configfile ../config_qc.yaml full/f1.pgen -j 10
```

:::
:::{tab-item} Sandbox

```bash
cd GDCGenomicsQC/workflow
gdcgenomicsqc --configfile ../config_qc.yaml full/f1.pgen -j 10
```

:::
:::{tab-item} Other HPCs

```bash
cd GDCGenomicsQC/workflow
gdcgenomicsqc --configfile ../config_qc.yaml full/f1.pgen -j 10
```

:::
:::{tab-item} Local Snakemake

```bash
cd GDCGenomicsQC/workflow
snakemake --profile=../profiles/hpc \
    --configfile ../config_qc.yaml \
    full/f1.pgen \
    -j 10
```

:::
::::

The Initial QC stage performs:

1. **Sample missingness (initial)**: Removes samples with >10% missing genotypes (``--mind 0.1``)
2. **Variant missingness**: Removes variants with >2% missingness (``--geno 0.02``)
3. **Sample missingness (final)**: Removes samples with >2% missingness (``--mind 0.02``)
4. **LD pruning**: Creates pruned dataset for downstream analyses (``--indep-pairwise 500 10 0.1``)

### Step 3: Run Standard QC

::::{tab-set}
:::{tab-item} MSI HPC

```bash
gdcgenomicsqc --configfile ../config_qc.yaml full/f1.b38.f2.pgen -j 10
```

:::
:::{tab-item} Sandbox

```bash
gdcgenomicsqc --configfile ../config_qc.yaml full/f1.b38.f2.pgen -j 10
```

:::
:::{tab-item} Other HPCs

```bash
gdcgenomicsqc --configfile ../config_qc.yaml full/f1.b38.f2.pgen -j 10
```

:::
:::{tab-item} Local Snakemake

```bash
snakemake --profile=../profiles/hpc \
    --configfile ../config_qc.yaml \
    full/f1.b38.f2.pgen \
    -j 10
```

:::
::::

The Standard QC stage applies additional filters:

1. **Minor Allele Frequency (MAF)**: Removes variants with MAF < 1% (``--maf 0.01``)
2. **Hardy-Weinberg Equilibrium (HWE)**: Removes variants failing HWE at p < 1×10⁻⁶ (discovery) and p < 1×10⁻¹⁰ (validation)
3. **Heterozygosity check**: Identifies samples with FWER > 3 standard deviations from mean
4. **Sex check**: Optionally verifies reported sex matches genetic sex

### Step 4: Run Ancestry-Specific QC

After ancestry classification, run QC on specific ancestry groups:

::::{tab-set}
:::{tab-item} MSI HPC

```bash
gdcgenomicsqc --configfile ../config_qc.yaml EUR/f1.b38.f2.pgen -j 10
```

:::
:::{tab-item} Sandbox

```bash
gdcgenomicsqc --configfile ../config_qc.yaml EUR/f1.b38.f2.pgen -j 10
```

:::
:::{tab-item} Other HPCs

```bash
gdcgenomicsqc --configfile ../config_qc.yaml EUR/f1.b38.f2.pgen -j 10
```

:::
:::{tab-item} Local Snakemake

```bash
snakemake --profile=../profiles/hpc \
    --configfile ../config_qc.yaml \
    EUR/f1.b38.f2.pgen \
    -j 10
```

:::
::::

Available subsets are dynamically determined from classification results.

### Visualizations

**Reference Space (PCA)**: ``images/PC_reference_space.svg``

- Reference panel samples in PC space with density contours
- Shows how target samples map to known ancestry groups

**Reference Space (UMAP)**: ``images/UMAP_reference_space.svg``

- Reference panel samples in UMAP embedding with density contours
- Nonlinear visualization of ancestry structure

---

## Interpreting Pipeline Outputs

### Sample Missingness Plot

**File**: ``{subset}/figures/smiss.svg``

![](images/smiss.svg)

The sample missingness histogram shows the distribution of missing data per individual.

- X-axis: Percentage of genotype calls missing per sample
- Red vertical lines: Threshold cutoffs (10% initial, 2% final)
- Samples to the right of the rightmost line are removed

**Note**: Since we used synthetic and (in case of the R25 data) imputed data, we don't expect to see any missingness in this exercise.

Standard interpretation:

- **Sharp peak at low values**: Good quality data
- **Long right tail**: Problematic samples requiring investigation
- **Bimodal distribution**: Possible batch effects or technology issues

### Variant Missingness Plot

**File**: ``{subset}/figures/vmiss.svg``

![](images/vmiss.svg)

The variant missingness histogram shows the distribution of missing data per SNP.

- X-axis: Percentage of samples missing genotype call per variant
- Red vertical lines: Threshold cutoffs
- Variants to the right of the line are removed

**Note**: Since we used synthetic and (in case of the R25 data) imputed data, we don't expect to see any missingness in this exercise.

Standard interpretation:

- **Concentrated at low values**: Clean variant calling
- **Long right tail**: Possible strand flip issues, poor quality regions, or structural variants

### Unplotted Output Files

The QC pipeline generates many intermediate files for detailed analysis:

**Initial QC Outputs**:

```{list-table}
:header-rows: 1

* - File
  - Description
* - ``f1.pgen/.pvar/.psam``
  - Merged, filtered dataset
* - ``f1.ldpruned.*``
  - LD-pruned dataset for PCA/relatedness
* - ``initial.smiss``
  - Sample missingness table
* - ``initial.vmiss``
  - Variant missingness table
```

**Standard QC Outputs**:

```{list-table}
:header-rows: 1

* - File
  - Description
* - ``f1.b38.f2.pgen/.pvar/.psam``
  - Final filtered dataset
* - ``f1.b38.f2.ldpruned.*``
  - Final LD-pruned dataset
* - ``MAF_check.afreq``
  - Allele frequency table
* - ``zoomhwe.hwe``
  - Variants failing HWE p < 1×10⁻⁵
* - ``indepSNP.prune.in``
  - Independent SNPs for heterozygosity
* - ``R_check.het``
  - Heterozygosity rate per sample
* - ``fail-het-qc.txt``
  - Samples failing heterozygosity filter
* - ``sex_discrepancy.txt``
  - Samples with sex check problems
```

**Sample contents of** ``initial.smiss``:

```{list-table}
:header-rows: 1

* - IID
  - FID
  - F_MISS
  - N_MISS
* - S001
  - FAM001
  - 0.001
  - 150
* - S002
  - FAM001
  - 0.008
  - 1200
```

**Sample contents of** ``R_check.het``:

```{list-table}
:header-rows: 1

* - IID
  - FID
  - O.HOM.
  - N.NM.
* - S001
  - FAM001
  - 2500
  - 3000
* - S002
  - FAM001
  - 2450
  - 3000
```

The heterozygosity rate is calculated as: ``(N.NM. - O.HOM.) / N.NM.``

---

## Discussion Points: Multi-Ancestry and Admixed Study Designs

These questions explore QC considerations for diverse and admixed populations:

1. **Ancestry-specific allele frequencies**: The MAF filter (default 1%) may remove
   informative variants in population-specific contexts. How should MAF thresholds
   differ between ancestry groups? Should multi-ancestry studies use uniform or
   group-specific thresholds?

2. **HWE assumptions**: HWE testing assumes a randomly mating population. For
   admixed individuals, this assumption is violated. Should HWE filters be applied
   before or after ancestry classification? How do systematic departures from HWE
   in admixed populations affect downstream analysis?

3. **Heterozygosity in admixed samples**: Admixed individuals have higher
   heterozygosity than homogeneous populations. Does the 3-SD threshold appropriately
   capture excess heterozygosity as an outlier versus natural admixture? How does
   this affect false positive rates?

4. **Missingness patterns**: Samples with high global ancestry proportions may
   have higher missingness if reference panels poorly represent their ancestry.
   How should missingness thresholds account for reference panel coverage across
   diverse populations?

5. **Sex chromosome handling**: The pseudoautosomal regions and sex chromosome
   ploidy differ between ancestries. How should X-chromosome heterozygosity filters
   be adjusted for multi-ancestry studies?

6. **Relatedness in family-structured populations**: For studies with family
   structure across ancestries, should KING or PRIMUS be used? How does population
   structure affect kinship coefficient estimates?

7. **Differential QC power**: Some ancestry groups may have more variants removed
   due to technology bias (e.g., array density). How does differential QC success
   affect downstream GWAS power and potential for bias?

8. **Strand alignment**: Poorly aligned variants show as missing in specific
   ancestry groups. How do you distinguish true missingness from strand issues
   in multi-ancestry data?

For theoretical foundations—including population genetics principles, statistical
tests for QC metrics, and best practices for diverse populations—refer to the
accompanying lecture materials.

---

## Data Types Reference

This section provides a high-level overview of the file formats used throughout
the QC pipeline. Understanding these formats helps interpret pipeline outputs
and troubleshoot data issues.

### Genotype Data Formats

**VCF/VCF.gz (Variant Call Format)**

The standard format for storing genetic variants. Files ending in ``.vcf.gz``
are compressed using bgzip for efficiency. VCF files contain:

- Header lines describing the format and reference genome
- Metadata lines describing FILTER, INFO, and FORMAT fields
- Data lines with chromosome, position, ID, reference allele, alternate allele,
  quality, filter status, and INFO fields
- Genotype calls in the sample columns

Each VCF file is typically indexed by a corresponding ``.csi`` file (see below)
to enable random access to specific genomic regions.

**PLINK BED/BIM/FAM**

The traditional PLINK binary genotype format. A complete dataset consists of:

- ``.bed`` file: Binary genotype matrix (major mode)
- ``.bim`` file: Variant information (chromosome, rs ID, genetic position,
  base pair position, allele 1, allele 2)
- ``.fam`` file: Family/sample information (family ID, individual ID, paternal
  ID, maternal ID, sex, phenotype)

The ``.bim`` and ``.fam`` files are plain text; only the genotype matrix is
stored in binary format for compactness.

**PLINK PGEN**

The newer PLINK2 binary genotype format, offering advantages over BED:

- ``.pgen`` file: Binary genotype matrix with flexible encoding
- ``.pvar`` file: Variant information (replaces ``.bim``)
- ``.psam`` file: Sample information (replaces ``.fam``)

PGEN supports multiple variant encoding schemes within a single file, enabling
more efficient storage and faster operations on large datasets. The pipeline
uses PGEN internally for all computations.

### Genomic Interval Files

**BED (Genomic Intervals)**

Browser Extensible Display format for genomic intervals. Standard 3-column BED:

- Column 1: Chromosome
- Column 2: Start position (0-based)
- Column 3: End position (1-based)

Extended BED formats (6-column, 12-column) include additional fields such
as name, score, strand, thick start/end, item RGB, and block counts/sizes.
Used for genomic region definitions, target files, and annotation tracks.

### Index Files

**CSI Index Files**

Tabix-style index files (``.csi``) that enable random access to compressed
VCF files. CSI (Coordinate Sorted Index) is an improvement over the older
TBI format, supporting files with many more variants (up to ~4 billion).

When you see ``chr22.vcf.gz.csi``, this index file allows tools to quickly
locate variants in a specific genomic region without scanning the entire file.

### Statistical and Output Files

**GRM (Genetic Relationship Matrix)**

Binary matrix storing pairwise genetic relationships between samples. Used
for estimating heritability and controlling for population structure. The
GRM is computed from QC-filtered variants and stored in PLINK2's ``.grm``
format with accompanying ``.grm.id`` sample mapping.

**KING Output Files**

Relatedness estimates computed by the KING software. Output includes:

- ``.king`` file: Kinship coefficients between sample pairs
- ``.king.id``: Sample identifiers corresponding to matrix rows/columns

Kinship coefficients are used to identify and remove related individuals,
with a default cutoff of 0.0884 (equivalent to 3rd-degree relatives).

**IDS and ID Mapping Files**

Plain text files listing sample or variant identifiers:

- ``.ids`` files: Lists of sample IDs for subsetting operations
- ``.psam``/``.fam`` files: Sample identifiers with optional family structure
- Used for sample selection, excluded sample lists, and cross-referencing
  between analysis stages

**Eigenvec/Eigenval Files**

Principal component analysis (PCA) results:

- ``.eigenvec``: Sample coordinates in PCA space (sample ID + PC values)
- ``.eigenval``: Variance explained by each principal component

Used for ancestry visualization, population structure correction, and
outlier detection. The pipeline generates these files for both internal
PCA and reference panel projection.

### Tabular Data Files

**TSV (Tab-Separated Values)**

Plain text format for structured data where columns are separated by tabs.
Pipeline outputs in TSV format include:

- ``.smiss``: Sample missingness statistics
- ``.vmiss``: Variant missingness statistics
- ``.afreq``: Allele frequency summaries
- ``.hwe``: Hardy-Weinberg equilibrium test results

TSV is preferred for data with variable-length fields or when compatibility
with Unix tools (cut, awk, sort) is desired.

**CSV (Comma-Separated Values)**

Plain text format for structured data with comma delimiters. Used for:

- Configuration files (YAML is preferred for pipeline config, but CSV appears
  in some auxiliary data files)
- Population labels and annotation files

CSV is simpler than TSV but can cause issues with fields containing commas;
TSV is generally preferred for genetic data.

**Plink Report Files**

Various PLINK output files in text format:

- ``.het``: Heterozygosity statistics per sample
- ``.hh``: Half-homozygous variant calls (potential Mendel errors)
- ``.irem``: Samples removed by quality filters

These files are plain text TSV format and can be inspected directly.

### File Format Compatibility

The pipeline automatically converts between formats as needed:

```{list-table}
:header-rows: 1

* - Input Format
  - Internal Format
  - Output Options
* - VCF.gz
  - PGEN
  - PGEN, VCF.gz
* - BED/BIM/FAM
  - PGEN
  - PGEN, BED
* - PGEN
  - PGEN
  - PGEN, VCF.gz
```

All intermediate computations use PGEN for efficiency. Final outputs can be
converted to VCF for compatibility with downstream tools.

---

## Next Steps

After completing this tutorial, proceed to:

- [](tutorial_ancestry_classification.md) - Classify ancestry using the QC-filtered data
- [](tutorial_heritability.md) - Estimate heritability using QC-filtered genotypes

## **Lab Materials**

- [QC Visualization Lab (Quarto)](labs/lab02_qc_visualization.qmd) - Interactive R notebook for visualizing QC outputs

**See also:**

- [](installation.md) - Software setup (if not already done)
- [](usage.md) - Running the full pipeline
- [](genomics.md) - Technical details on QC methods
