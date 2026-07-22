# Genomics

This section outlines the standard procedures for the genomic data processing pipeline
and how to configure them via the config file.

## Standard Procedure

### Module 1: Crossmap (optional)
CrossMap converts genome coordinates and annotation files between different reference assemblies.
It supports a wide array of file formats, including BAM, CRAM, SAM, VCF, Wiggle, BigWig,
BED, GFF, and GTF. For the purposes of our pipeline we will make use of PLINK Binary
format and convert the genome build to GRCh38 (default: GRCh37 to GRCh38).

- **Config option:** `liftover: true` (Default: `true`)
- **Description:** Set to `false` if the input data is already GRCh38.

Example in `config.yaml`:

```yaml
liftover: true  # Set to false if data is already GRCh38
```

### Module 2: GenotypeHarmonizer (optional)
GenotypeHarmonizer integrates genetic data by resolving inconsistencies in genomic strand
and file format. It will align study datasets to a specified reference genome and uses
linkage disequilibrium (LD) to solve unknown or ambiguous strand issues and SNPs. We
use it in our pipeline to align sample data to the GRCh38 reference genome using the
reference file `ALL.hgdp1kg.filtered.SNV_INDEL.38.phased.shapeit5.vcf`.

- **Config option:** `harmonize: true` (Default: `true`)
- **Description:** Set to `false` if the data is already harmonized.

Example in `config.yaml`:

```yaml
harmonize: true  # Set to false if data is already harmonized
```

### Module 3: Initial QC
Performs Quality Control prior to relatedness checks.

- Exclude SNPs and individuals with >10% missingness (**Plink**; configurable via ``initial_variant_missingness`` and ``initial_subject_missingness``).
- Exclude SNPs and individuals with >2% missingness (**Plink**; configurable via ``final_variant_missingness`` and ``final_subject_missingness``).

The variant missingness filter is applied in two stages:

1. ``initial_variant_missingness`` (default: 0.1) — applied per-chromosome in the ``convertPlink`` step, before allele alignment.
2. ``final_variant_missingness`` (default: 0.02) — applied in ``initialFilter.sh`` after allele alignment, using a more stringent threshold.

Sample missingness is also filtered in two stages within ``initialFilter.sh``:

1. ``initial_subject_missingness`` (default: 0.1) — applied before the final variant missingness filter.
2. ``final_subject_missingness`` (default: 0.02) — applied after the final variant missingness filter.

All four thresholds can be overridden in ``config.yaml`` and are ``nullable`` (set to ``null`` to skip the filter).

### Module 4: Relatedness
This module uses **KING** to separate related and unrelated study samples. The module
also performs **PC-AiR** and **PC-Relate** for kinship estimation.

**KING** is a toolset that makes use of SNP data to identify how closely two individuals
are related based on their DNA. This inference is based on the **kinship coefficient** ($\phi$).

The fundamental equation KING uses to estimate the kinship coefficient $\phi$ between
individuals $i$ and $j$ is based on the counts of **Heterozygote-Heterozygote (Het-Het)**
and **Heterozygote-Homozygote (Het-Hom)** mismatches:

$$
\phi_{ij} = \frac{N_{Aa,Aa} - 2N_{AA,aa}}{N_{Aa,i} + N_{Aa,j}}
$$

Based on the calculated $\phi$, relationships are categorized as follows:

```{list-table} Relationship Inference Criteria
:widths: 30 30 40
:header-rows: 1

* - Relationship Degree
  - Kinship Coefficient ($\phi$)
  - Examples
* - **Duplicate/Twin**
  - $> 0.354$
  - Identical twins, same person sampled
* - **1st-Degree**
  - $[0.177, 0.354]$
  - Parent-Offspring, Full Siblings
* - **2nd-Degree**
  - $[0.0884, 0.177]$
  - Grandparent-Grandchild, Half-siblings
* - **Unrelated**
  - $< 0.0442$
  - No close detectable relation
```

PC-AiR is used to perform Principal Component Analysis (PCA) for population structure
detection while accounting for known or cryptic relatedness. It identifies a subset of
unrelated individuals that represent the ancestral diversity of the sample to compute
the principal components (PCs).

The method utilizes the kinship coefficients ($\phi$) calculated by **KING** to
define a "partition" of the data.

PC-Relate uses the principal components from PC-AiR to estimate kinship coefficients
and IBS (Identity By State) sharing probabilities while adjusting for population stratification.

The primary metric produced is the **PC-Relate Kinship Coefficient** ($\phi_{ij}^{PCR}$),
which is estimated using a ratio of genetic covariance adjusted for local ancestry:

$$
\phi_{ij}^{PCR} = \frac{\text{Cov}(G_i, G_j)}{2 \sqrt{\text{Var}(G_i) \text{Var}(G_j)}}
$$

Config Options:

```yaml
relatedness:
    method: "king"  # Options: "king", "primus", or other (assumes unrelated)
    king_cutoff: 0.0884  # KING cutoff for "unrelated"

SEX_CHECK: true  # Whether to perform sex check
GRM: true  # Whether to compute GRM
```

```{list-table} Relatedness Method Comparison
:widths: 20 40 40
:header-rows: 1

* - Method
  - Primary Strength
  - Usage in Pipeline
* - **KING**
  - Robust to population structure without needing PCs.
  - Initial relatedness screening and PC-AiR partitioning.
* - **PC-AiR**
  - Captures ancestry without bias from family clusters.
  - Generating ancestry PCs for regression models.
* - **PC-Relate**
  - High accuracy in admixed populations.
  - Final kinship estimation and relatedness filtering.
```

## Internal PCA Methods

The pipeline supports two methods for computing internal PCA:

```{list-table} Internal PCA Method Comparison
:widths: 20 40 40
:header-rows: 1

* - Method
  - Description
  - Output
* - **plink2**
  - Fast approximate PCA using PLINK2 on unrelated samples
  - `internal_pca_plink2.eigenvec`, `internal_pca_plink2.eigenval`
* - **pcair**
  - PC-AiR using all samples, computes PC-relate kinship and GRM
  - `pcair_pcaobj.RDS`, `pcrelate_kinship.RDS`, `pcair.grm.bin`
* - **both**
  - Run both methods
  - `internal_pca_plink2.eigenvec`, `pcair_pcaobj.RDS`
```

```yaml
internalPCA:
    method: "plink2"  # Options: "plink2", "pcair", "both"
    npc: 20           # Number of PCs for plink2 method
```

### Module 5: Standard QC
Standard GWAS quality control measures on unrelated individuals.

- **Filtering:** Exclude SNPs and individuals with >2% missingness (**Plink**).
- **MAF:** Exclude SNPs with Minor Allele Frequency < 0.01.
- **HWE:** Exclude SNPs with p-values < 1e-6 (controls) or < 1e-10 (cases). The threshold can be scaled by sample size via the ``hwe_k`` config parameter (default: null, meaning fixed threshold). When set, the effective p-value is p × 10^(−n×k). Greer et al. (2024) recommends k=0.001 for large studies.
- **Sex Check:** F-values < 0.2 assigned as female, > 0.8 as male.

### Module 6: Phasing
Phasing performed via **shapeit4.2** with reference map `chr${CHR}.b38.gmap.gz`.

Phasing is the process of estimating haplotypes from observed genotypes, determining
which alleles were inherited together from a single parent. In this pipeline, phasing
is performed via **shapeit4.2**.

The accuracy of the haplotype estimation relies on a high-resolution genetic map that
provides the recombination rates across the genome. We use the reference map: `chr${CHR}.b38.gmap.gz`.

A common metric for evaluating phasing quality is the **Switch Error Rate**, which
measures the frequency of incorrect "switches" between the maternal and paternal
haplotypes in the estimated sequence:

$$
\text{SER} = \frac{\text{Number of Switch Errors}}{\text{Total Number of Opportunities for Switch Errors}}
$$

```{list-table} Phasing Parameters and Resources
:widths: 30 70
:header-rows: 1

* - Parameter/Resource
  - Description
* - **Software**
  - **shapeit4.2**: A fast and accurate method for estimation of haplotypes.
* - **Reference Map**
  - `chr${CHR}.b38.gmap.gz`: Genetic map used to model recombination.
* - **Input Format**
  - VCF/BCF: Requires high-quality, QC-filtered genotypes from Module 6.
* - **Output Format**
  - Phased VCF: Necessary for local ancestry inference in Module 7.
```

### Module 7: Rfmix
This module infers local ancestry across the genome using phased genotype files and a
reference panel, such as `hg38_phased.vcf.gz`. **rfmix** uses a discriminative
machine learning approach to assign ancestral origins to specific chromosomal segments.

For high-confidence ancestry calls, the pipeline enforces a strict threshold on the
posterior probabilities assigned to each segment. Global ancestry estimates are only
calculated for individuals where the posterior probability exceeds 0.8.

The posterior probability $P(A | G)$ represents the likelihood that a genomic
segment belongs to ancestry $A$ given the observed genotypes $G$:

$$
P(A | G) = \frac{P(G | A) P(A)}{P(G)}
$$

Config Options:

```yaml
localAncestry:
    RFMIX: true  # Enable RFMix
    test: true   # Run in test mode with reduced parameters
    thin_subjects: 0.1  # Fraction of subjects to use
```

```{list-table} rfmix Configuration and Requirements
:widths: 30 70
:header-rows: 1

* - Parameter/Resource
  - Requirement/Description
* - **Input Files**
  - Must be phased VCF files from **shapeit4.2** (Module 6).
* - **Reference Panel**
  - `hg38_phased.vcf.gz`: Phased reference genotypes for known populations.
* - **Genetic Map**
  - Requires a genetic map (recombination rates) consistent with the genome build.
* - **Confidence Threshold**
  - Posterior probability $> 0.8$ for global ancestry inclusion.
```

### Module 8: Global Ancestry Classification
This module classifies global ancestry using dimension reduction (PCA, UMAP, VAE)
and Random Forest classification. See [](tutorial_ancestry_classification.md) for detailed usage.

Config Options:

```yaml
ancestry:
    threshold: 0.8  # Minimum posterior probability for confident classification
    model: "pca"    # Embedding model: pca, umap, vae, rfmix
```

## File Naming Conventions

The pipeline uses a structured naming convention for intermediate and output files.
Each stage appends a suffix to the file prefix indicating what processing steps have
been applied:

```{list-table} Naming Convention Suffixes
:widths: 25 75
:header-rows: 1

* - Suffix
  - Meaning
* - `f1`
  - **Filter 1**: Initial QC — sample and variant missingness filtering (`--geno`, `--mind`)
* - `f1.f2`
  - **Filter 2**: Standard QC — MAF, HWE, heterozygosity, and optional sex check applied on top of f1
* - `f1.ldpruned`
  - **LD Pruned**: Variants pruned for linkage disequilibrium using `--indep-pairwise 50 5 0.2`
* - `f1.ldpruned.unrelated`
  - **Unrelated subset**: Related samples removed via KING (or other method) from the LD-pruned set
* - `f1.ldpruned.unrelated.ldpruned`
  - **LD Pruned (unrelated)**: A second LD-pruning pass applied to the unrelated subset
* - `f1.f2.ldpruned`
  - **LD Pruned (post-Standard QC)**: LD pruning applied after both QC filters
```

Example output file tree for the `full` subset:

```
OUT_DIR/
└── full/
    ├── f1.pgen                        # After initial QC
    ├── f1.f2.pgen                     # After standard QC
    ├── f1.ldpruned.pgen               # LD pruned (pre-relatedness)
    ├── f1.ldpruned.unrelated.pgen     # Unrelated samples only
    ├── f1.ldpruned.unrelated.ldpruned.pgen  # LD pruned unrelated samples
    ├── f1.f2.ldpruned.pgen            # LD pruned post-standard QC
    ├── f1.ldpruned.unrelated.grm.bin  # GRM from unrelated subset
    ├── f1.ldpruned.unrelated.grm.id
    ├── f1.ldpruned.unrelated.grm.N.bin
    ├── internal_pca_plink2.eigenvec   # Internal PCA (plink2)
    └── internal_pca_plink2.eigenval
```

When processing per-chromosome inputs, chromosome-specific files use the pattern
`f1_{CHR}.pgen` and `f1.f2_{CHR}.pgen`.

## Technical Implementation

### Module 1: Crossmap
```console
$ python CrossMap.py GRCh37_to_GRCh38.chain.gz prep.bed study_lifted.bed
```

### Module 2: GenotypeHarmonizer
```console
$ java -jar GenotypeHarmonizer.jar --input study_data --input_type PLINK_BED \
  --ref ALL.hgdp1kg.filtered.SNV_INDEL.38.phased.shapeit5.vcf --ref_type VCF \
  --output harmonized_data --output_type PLINK_BED
```

### Module 3: Initial QC
```console
$ plink --bfile file_stem --geno 0.02 --make-bed QC1
$ plink --bfile QC1 --mind 0.02 --make-bed --out QC2
```

### Module 4: KING and PC-AiR
```console
# Run KING to get kinship estimates
$ king -b study.bed --kinship --prefix king_results

# Example R code snippet for PC-AiR/PC-Relate via GENESIS
$ Rscript run_genesis.R --king king_results.kin0 --vcf study.vcf.gz
```

### Module 5: Standard QC (Sex Check)
```console
$ plink --bfile QC4 --check-sex
$ grep 'PROBLEM' plink.sexcheck | awk '{print $1, $2}' > sex_discrepancy.txt
```

### Module 6: Phasing
```console
# Execute shapeit4 for a specific chromosome
$ shapeit4 --input study_filtered.vcf.gz \
           --map chr${CHR}.b38.gmap.gz \
           --region ${CHR} \
           --output study_phased_chr${CHR}.vcf.gz \
           --thread 8
```

### Module 7: Rfmix Execution
```console
# Execute rfmix for local ancestry inference
$ rfmix -f study.phased.vcf.gz \
        -r reference_panel.phased.vcf.gz \
        -m sample_map.txt \
        -g genetic_map.txt \
        -e 2 \
        -n 5 \
        --chromosome=chr${CHR} \
        -o output_ancestry
```
