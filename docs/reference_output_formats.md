---
title: "Pipeline Output File Reference"
---

# Overview

This page catalogs every distinct output format produced by the GDCGenomicsQC
pipeline, organized by file extension. Use it to find premade tools for standard
operations and the R/Python code needed to read each file type for methods
research.

Abbreviations used in the table:

- **Ext** — file extension(s)
- **Info** — what the file contains
- **Use** — common downstream research applications
- **Tools** — pre-built software that operates on this format
- **R** — how to load into R
- **Py** — how to load into Python

# Genome formats

| Ext | Info | Use | Tools | R | Py |
|-----|------|-----|-------|----|----|
| `.pgen` | Binary genotype matrix | GWAS, PRS, phasing, PCA, relatedness | PLINK2, PRS-CSx, LDpred2, GCTA (`--mgrm`), BOLT-LMM, REGENIE, SAIGE | `genio::read_plink("prefix")` returns a list with `$X` (genotypes as `FBM` via `bigsnpr`), `$fam`, `$map`; `bigsnpr::snp_readBed2()` for PLINK1; `gaston::read.bed.matrix()` | `pgenlib.Pgen` (Python `pgenlib` crate); `bed_reader.open_bed()` reads PLINK1 `.bed`; `pandas_plink.read_plink()` loads into xarray |
| `.pvar` | Tabular variant metadata. Columns: `#CHROM` `POS` `ID` `REF` `ALT` | Variant lookup, filtering by gene/region, MAF extraction, tallying by chromosome | PLINK2, BCFtools/py (intersect with VCF regions), `rtracklayer` (import as GRanges) | `readr::read_tsv("file.pvar")`; `data.table::fread()`; `rtracklayer::import()` as GRanges | `pandas.read_table("file.pvar")`; `polars.read_csv("file.pvar", separator="\t")` |
| `.psam` | Tabular sample metadata. Columns: `#IID` `FID` `SEX` `PHENO` | Sample lookup, phenotype extraction, sex filter | PLINK2, any text tool | `readr::read_tsv("file.psam")`; `data.table::fread()` | `pandas.read_table("file.psam")` |
|-----|------|-----|-------|----|----|
| `.bed` (PLINK1) | Binary genotype matrix | GWAS, PRS, clumping, LD estimation, PCA. Most legacy tools still require PLINK1 format. | PLINK1.9, PRS-CS, LDpred, GCTA, EIGENSOFT, ADMIXTURE, KING, BOLT-LMM, REGENIE, SAIGE | `genio::read_plink("prefix")`; `bigsnpr::snp_readBed2("file.bed", backingfile="tmp")` → `FBM`; `gaston::read.bed.matrix("prefix")`; `SNPRelate::snpgdsBED2GDS()` | `bed_reader.open_bed("file.bed")` → numpy array; `pandas_plink.read_plink("prefix")`; `scikit-allel.read_vcf()` for VCF-derived |
| `.bim` | Extended variant map. Columns: `CHR` `SNP` `GD` `BP` `A1` `A2` | LD-pruning, clumping, variant QC, MAF calculation, chromosome filtering | PLINK1.9, any text tool | `readr::read_tsv("file.bim", col_names=F)`; `genio::read_plink("prefix")$map` | `pandas.read_table("file.bim", header=None)` |
| `.fam` | Sample pedigree + phenotype. Columns: `FID` `IID` `PAT` `MAT` `SEX` `PHENO` | Sample filtering, phenotype extraction, case/control designation | PLINK1.9, any text tool; `gap::GWAS.by.chr()` | `readr::read_tsv("file.fam", col_names=F)`; `genio::read_plink("prefix")$fam` | `pandas.read_table("file.fam", header=None)` |
|-----|------|-----|-------|----|----|
| `.vcf.gz` | BGZF-compressed VCF with genotype likelihoods or phased calls. INFO, FORMAT, GT fields vary by stage. | Phasing input/output, imputation (IMPUTE5/minimac4), local ancestry (RFMix), population genetics (FST, Tajima's D), GWAS (dosage-based), fine-mapping | BCFtools, Shapeit5, IMPUTE5, minimac4, Beagle5, whatsHap, vcftools, PLINK2 (`--vcf`), bcftools, htslib, GATK | `SeqArray::seqVCF2GDS("in.vcf.gz", "out.gds")` then `SeqArray::seqOpen()`; `vcfR::read.vcf("file.vcf.gz")` → `vcfR` object; `gwasvcf::readVcf()`; `snpStats::read.pedfile()` for PLINK text | `pysam.VariantFile("file.vcf.gz")` → iterable; `cyvcf2.VCF("file.vcf.gz")`; `scikit-allel.read_vcf("file.vcf.gz")`; `sgkit.io.vcf.read_vcf()` |
| `.vcf.gz.csi` | CSI index for BGZF-compressed VCF/BCF. Enables efficient random access by genomic region. | Efficient region-based queries (e.g., extract chr22:1-50000) | BCFtools, htslib, pysam | `pysam.VariantFile("file.vcf.gz")` auto-loads `.csi` if present | `pysam.VariantFile("file.vcf.gz")` auto-index |
| `.vcf.gz.tbi` | Tabix index (alternative to `.csi`). Same purpose, different format. | Region lookups, IGV loading | tabix (htslib), IGV, BCFtools | `pysam.VariantFile("file.vcf.gz")` also detects `.tbi` | `pysam.VariantFile()` auto-detects `.tbi`; `cyvcf2.VCF()` |



## Relatedness formats

| Ext | Info | Use | Tools | R | Py |
|-----|------|-----|-------|----|----|
| `.grm.bin` | Lower-triangular GRM as flat binary (float32, N×N / 2 entries). N = # samples from `.grm.id`. | SNP heritability (GREML), genetic correlation, PCA (Eigenstrat-style), variance decomposition, MTAG | GCTA, GCTB, BOLT-REML, BOLT-LMM, MTG2, LDAK, MEGHA, LMR | `N <- nrow(id); m <- N*(N+1)/2; g <- readBin("file.grm.bin", double(), m, size=4)` then reconstruct lower triangle | `import numpy as np; N = len(id); g = np.fromfile("file.grm.bin", dtype=np.float32);` reconstruct lower triangle via `np.tril_indices` |
| `.grm.id` | Sample order of GRM (FID+IID, one per line, N lines). | Maps GRM rows/cols to sample IDs | Any text tool | `readr::read_table("file.grm.id", col_names=F)`; `data.table::fread()` | `pandas.read_table("file.grm.id", header=None)` |
| `.grm.N.bin` | Per-pair non-missing SNP counts (float64, same layout as `.grm.bin`). | QC — identifies poorly-estimated GRM entries (low N), filtering | GCTA | `readBin("file.grm.N.bin", double(), m)` | `np.fromfile("file.grm.N.bin", dtype=np.float64)` |
|-----|------|-----|-------|----|----|
| `.kin0` | KING pairwise kinship coefficients. Columns: `FID1` `IID1` `FID2` `IID2` `N_SNP` `Kinship` | Relatedness inference, sample exclusion (KING cutoff ≥0.0884 = 3rd-degree relatives), pedigree validation | KING, PLINK2 (`--make-king`), `SNPRelate::snpgdsIBD()` | `readr::read_table("file.kin0")`; `dplyr::filter(kin, Kinship > 0.0884)`; `coxme::kinship()` | `pandas.read_table("file.kin0")` |
| `.genome` | PLINK1 IBD estimation (PI_HAT). Columns: `FID1` `IID1` `FID2` `IID2` `RT` `EZ` `Z0` `Z1` `Z2` `PI_HAT` `PHE` `DST` `PPC` `RATIO` | IBD-based relatedness (alternative to KING), PRIMUS input, duplicate/ MZ twin detection | PLINK1.9 (`--genome`), PRIMUS, KING | `readr::read_table("file.genome")`; `dplyr::filter(genome, PI_HAT > 0.125)` | `pandas.read_table("file.genome")` |
| `*.king.cutoff.in.id` | Sample IDs within KING kinship cutoff (unrelated set). One FID+IID per line. | Downstream analysis that requires unrelated samples (PCA, heritability, GWAS) | Any text tool; `--keep` with PLINK2 | `readr::read_table("file.king.cutoff.in.id", col_names=F)` | `pandas.read_table("file.king.cutoff.in.id", header=None)` |
| `*.king.cutoff.out.id` | Sample IDs exceeding KING kinship cutoff (related set). | Related samples list, pedigree QC | Any text tool | `readr::read_table("file.king.cutoff.out.id", col_names=F)` | `pandas.read_table("file.king.cutoff.out.id", header=None)` |
|-----|------|-----|-------|----|----|
| `.eigenvec` | PCA sample coordinates. Columns: `#IID` `FID` `PC1` `PC2` ... `PC{N}` | Population structure correction (covariates in GWAS), ancestry inference, outlier detection, visualisation | PLINK2, EIGENSOFT (`smartpca`), FlashPCA, GCTA, TeraPCA | `readr::read_table("file.eigenvec")`; `ggplot(aes(x=PC1, y=PC2))` | `pandas.read_table("file.eigenvec")`; `plt.scatter(df.PC1, df.PC2)` |
| `.eigenval` | PCA eigenvalues, one per line. Variance explained by each PC. | Determine number of significant PCs (scree plot, elbow test) | PLINK2, any text tool | `scan("file.eigenval")`; `plot(prop.table(val))` | `np.loadtxt("file.eigenval")`; `plt.plot(val/sum(val))` |
| `.eigenvec.allele` | PCA loading / SNP weights. Columns: `#CHROM` `POS` `ID` `REF` `ALT` `PC1_WT` `PC2_WT` ... | Identify variants driving PC axes (genomic annotation of population structure) | PLINK2 (`--pca allele-wts`), `pca_loading_plot()` in custom R | `readr::read_table("file.eigenvec.allele")`; `dplyr::top_n(loadings, 10, abs(PC1_WT))` | `pandas.read_table("file.eigenvec.allele")` |
| `.sscore` | PLINK2 polygenic score / projection output. Columns: `#IID` `FID` `SCORE1_AVG` `SCORE1_SUM` `NAMED_ALLELE_DOSAGE_SUM` (variable) | PRS calculation, PCA projection of study samples onto reference PCs | PLINK2 (`--score`, `--pca-proj`), any text tool | `readr::read_table("file.sscore")` | `pandas.read_table("file.sscore")` |

## Ancestry formats

| Ext | Info | Use | Tools | R | Py |
|-----|------|-----|-------|----|----|
| `classificationProbabilities.tsv` | Per-sample ancestry posterior probabilities. Columns: `IID`, one column per ancestry class (`AFR`, `EUR`, `EAS`, `AMR`, `SAS`, `Other`). | Ancestry-aware downstream analysis, admixture proportion estimation, population genetics | Any text tool; R: `heatmap(as.matrix(probs))`; Py: `sns.heatmap()` | `readr::read_tsv("classificationProbabilities.tsv")` | `pandas.read_table("classificationProbabilities.tsv")` |
| `ancestry_classifications.tsv` | Final ancestry labels per sample. Columns: `IID`, model predictions (`pca_predicted`, `umap_predicted`), `probability`. | Sample ancestry assignment, per-ancestry subsetting for stratified analysis | Any text tool; PLINK2 `--keep` with generated keep files | `readr::read_tsv("ancestry_classifications.tsv")`; `dplyr::filter(anc, pca_predicted=="EUR")` | `pandas.read_table("ancestry_classifications.tsv")` |
| `keep_{ANC}.txt` | One IID per line, one file per ancestry (e.g., `keep_EUR.txt`, `keep_AFR.txt`). | Sample filtering by ancestry; direct input to PLINK2 `--keep` | PLINK2 (`--keep`), any text tool | `readr::read_table("keep_EUR.txt", col_names="IID")` | `pandas.read_table("keep_EUR.txt", header=None, names=["IID"])` |
| `umap_ref.csv` / `umap_sample.csv` | UMAP 2D coordinates. Columns: `IID`, `UMAP1`, `UMAP2`. | Non-linear dimension reduction for ancestry visualization; complements PCA | Any text tool; R: `ggplot(aes(x=UMAP1, y=UMAP2))`; Py: `plt.scatter(df.UMAP1, df.UMAP2)` | `readr::read_csv("umap_ref.csv")` | `pandas.read_csv("umap_ref.csv")` |

## Local ancestry (RFMix) formats

| Ext | Info | Use | Tools | R | Py |
|-----|------|-----|-------|----|----|
| `*.lai.msp.tsv` | RFMix Most Recent Shared Ancestor calls. Columns: chromosome, start, end, then per-haplotype ancestry calls (2 columns per sample). | Locus-level ancestry assignment, ancestry-specific GWAS, fine-mapping with ancestry-aware LD, admixture mapping | RFMix, FLARE, LAMP-LD, `plot_rfmix.py` (RFMix utils) | `readr::read_tsv("file.lai.msp.tsv")`; reshape per-haplotype columns; `ggplot()` with `geom_rect()` for karyotype-style plots | `pandas.read_table("file.lai.msp.tsv")`; `matplotlib` or `plotnine` for visualization |
| `*.lai.fb.tsv` | RFMix forward-backward posterior probabilities. Columns: position, then per-ancestry probability per haplotype. | Confidence assessment of local ancestry calls, uncertainty-aware downstream analysis | RFMix, custom scripts | Read similarly to `.msp.tsv`; `rowSums()` to check calibration | Read similarly to `.msp.tsv`; `df.sum(axis=1)` to verify proper normalization |
| `*.lai.rfmix.Q` | Global ancestry proportions from RFMix (genome-wide average of local calls). First 2 header lines, then sample_id + per-ancestry proportions. | Genome-wide admixture proportions (alternative/validation of global ancestry classification) | Any text tool; compare to `classificationProbabilities.tsv` | `readr::read_tsv("file.lai.rfmix.Q", skip=2)` | `pandas.read_table("file.lai.rfmix.Q", skiprows=2)` |
| `ancestry_full.txt` | Genome-wide aggregated local ancestry proportions per sample across all chromosomes. | Final per-sample local ancestry proportions for downstream analysis | Any text tool | `readr::read_tsv("ancestry_full.txt")` | `pandas.read_table("ancestry_full.txt")` |

## Association testing formats

| Ext | Info | Use | Tools | R | Py |
|-----|------|-----|-------|----|----|
| `*.glm.linear` | PLINK2 linear regression results. Columns: `#CHROM` `POS` `ID` `REF` `ALT` `A1` `TEST` `OBS_CT` `BETA` `SE` `T_STAT` `P` | GWAS (quantitative trait), SNP-trait association, meta-analysis input | PLINK2, METAL, GWAMA, MR-MEGA, PLINK clump, LDSC (formatted), FUMA | `readr::read_tsv("file.glm.linear")`; `qqman::qq(gwas$P)`; `qqman::manhattan(gwas)`; `MAGMA` (reformat) | `pandas.read_table("file.glm.linear")`; `scikit-allel` manhattan; `plotly` interactive |
| `*.glm.logistic` | PLINK2 logistic regression results. Similar columns but with `OR` `CI_LOW` `CI_UP` `P` instead of `BETA` `SE`. | GWAS (binary/case-control trait) | Same as `.glm.linear` | `readr::read_tsv("file.glm.logistic")` | `pandas.read_table("file.glm.logistic")` |
| `*_sumstats.txt` / `*_sumstats_singlePRS.txt` | Formatted GWAS summary statistics. Columns vary by PRS method (`SNP` `CHR` `A1` `A2` `BETA` `SE` `P` `N` and `rsid` versions). | PRS calculation (PRS-CS, LDpred2, PRS-CSx, SDPR, lassosum2, CT-SLEB), meta-analysis | PRS-CS, PRS-CSx, LDpred2, SDPR, lassosum2, CT-SLEB, MTAG, LDSC, FINEMAP, susieR, COJO | `readr::read_table("file_sumstats.txt")`; `ieugwasr::format_data()`; `TwoSampleMR::read_exposure_data()` | `pandas.read_table("file_sumstats.txt")`; `scipy.stats` for derived statistics |

## Heritability formats

| Ext | Info | Use | Tools | R | Py |
|-----|------|-----|-------|----|----|
| `herit.csv` | SNP heritability estimates. Columns: `Ancestry` `Phenotype` `h2` `SE` `Method` `N` `PCs` (varies by method) | Heritability interpretation, power analysis, genetic architecture, GxE comparison | GCTA, BOLT-REML, MASH, LDSC, SUMHER, GREML | `readr::read_csv("herit.csv")`; `ggplot(aes(x=Ancestry, y=h2))` | `pandas.read_csv("herit.csv")`; `seaborn.barplot()` |

## Phenotype & simulation formats

| Ext | Info | Use | Tools | R | Py |
|-----|------|-----|-------|----|----|
| `.pheno` | Phenotype file. Columns: `FID` `IID` `PHENO` (with header). | GWAS phenotype input, case/control designation, quantitative trait analysis | PLINK2 (`--pheno`), any text tool | `readr::read_tsv("file.pheno")`; `gwasglue::gwas_data()` | `pandas.read_table("file.pheno")` |
| `simulated.fam` | PLINK1 `.fam` with simulated phenotypes in column 6. | Validate heritability estimation (compare input h2 to estimated h2), power analysis for PRS/GWAS | PLINK1.9, GCTA, any analysis software | `genio::read_plink("simulated")$fam`; `readr::read_table("simulated.fam", col_names=F)` | `pandas.read_table("simulated.fam", header=None)` |
| `simulated_pheno1.pheno` | Extracted phenotype column from simulation. Columns: `FID` `IID` `PHENO`. | Direct phenotype input for heritability/GWAS on simulated data | Same as `.pheno` | `readr::read_tsv("simulated_pheno1.pheno")` | `pandas.read_table("simulated_pheno1.pheno")` |

## PRS intermediate formats

| Ext | Info | Use | Tools | R | Py |
|-----|------|-----|-------|----|----|
| `{method}_manifest.tsv` | Key-value pairs of all input paths for a PRS method run. Columns: `KEY` `VALUE`. | Debugging PRS pipeline, understanding input file paths used by each PRS method | Any text tool | `readr::read_tsv("manifest.tsv")`; `tibble::deframe()` | `pandas.read_table("manifest.tsv")`; `df.set_index("KEY")["VALUE"]` |
| `prs_inputs.env` | KEY=VALUE environment variable file. Paths to all inputs required by the PRS pipelines. | Sourced by shell scripts; can be `source`'d in R or Python to access paths programmatically | Shell (`source`), any text tool | `read.table("prs_inputs.env", sep="=")` | `pandas.read_table("prs_inputs.env", sep="=", comment="#")` |
| `prs_prscsx_generated.conf` / `prs_single_ancestry_*_generated.conf` | KEY=VALUE configuration files for PRS-CSx and single-ancestry PRS methods. | PRS method configuration, reproducibility | PRS-CSx, PRS-CS, any text tool | `read.table("file.conf", sep="=") %>% deframe()` | `pandas.read_table("file.conf", sep="=", comment="#")` |

## R object formats

| Ext | Info | Use | Tools | R | Py |
|-----|------|-----|-------|----|----|
| `.RDS` | Arbitrary R serialized object (single). Can contain PCA model (`pcair_pcaobj.RDS`), kinship estimates (`pcrelate_kinship.RDS`), or RF classifier (`RFpc.Rds`). | Reload trained models without retraining, access PC-AiR/PC-Relate results for custom plotting | `base::readRDS()`, `base::saveRDS()`, `base::load()` (for `.RData`) | `readRDS("file.RDS")` | `pyreadr.read_r("file.RDS")` (limited, depends on object type); for PC-AiR: use `pcair_coordinates.tsv` instead |
| `.gds` / `.seq_gds` | SeqArray / CoreArray Genomic Data Structure. Binary, indexed, supports efficient read-write on disk. | Intermediate format for PC-AiR, PC-Relate; efficient variant lookup without loading entire genotype matrix | SNPRelate, SeqArray, GDSArray, GENESIS, ArrayTools | `SeqArray::seqOpen("file.gds")`; `SeqArray::seqGetData()`; `SNPRelate::snpgdsOpen()` | `gdsfmt` Python package (limited); use `SeqArray::seqGDS2VCF()` to convert then read with pysam |

## Other useful text formats

| Ext | Info | Use | Tools | R | Py |
|-----|------|-----|-------|----|----|
| `.prune.in` | LD-independent variant IDs (one per line). | Input to `--extract` for PCA, heritability, PRS; ensures markers are unlinked for GCTA GRM | PLINK2 (`--extract`), any text tool | `readr::read_table("file.prune.in", col_names="ID")` | `pandas.read_table("file.prune.in", header=None, names=["ID"])` |
| `.prune.out` | LD-pruned variant IDs (one per line). | Variants removed due to LD; can be re-examined for region-specific LD patterns | PLINK2, any text tool | `readr::read_table("file.prune.out", col_names="ID")` | `pandas.read_table("file.prune.out", header=None, names=["ID"])` |
| `flip_list.txt` | Variant IDs requiring strand flip (allele complement mismatch). | Allele alignment debugging, strand-aware PRS (ensure effect alleles match) | PLINK2 (`--flip`), any text tool | `readr::read_table("flip_list.txt", col_names="ID")` | `pandas.read_table("flip_list.txt", header=None, names=["ID"])` |
| `align_report.txt` | Allele alignment summary counts. Lines: `exact_match: N`, `strand_flip: N`, `allele_mismatch: N`, `study_only: N`, `ref_only: N`. | QC — check alignment quality between study data and reference | Any text tool | `read.delim("align_report.txt", sep=":", header=F)` | `pandas.read_csv("align_report.txt", sep=":", header=None)` |
| `zoomhwe.hwe` / `zoomhwe_{CHR}.hwe` | Variants with HWE p < 1e-5 (from awk filter). Same columns as `.hardy`. | Flag variants with suggestive HWE deviation for review | Any text tool; grep for specific variants | `readr::read_table("zoomhwe.hwe")` | `pandas.read_table("zoomhwe.hwe")` |
| `sex_discrepancy.txt` / `sex_discrepancy_{CHR}.txt` | FID IID of samples with mismatched sex. | Remove sex-mismatched samples before downstream analysis | PLINK2 (`--remove`), any text tool | `readr::read_table("sex_discrepancy.txt", col_names=F)` | `pandas.read_table("sex_discrepancy.txt", header=None)` |
| `het_fail_ind.txt` / `het_fail_clean.txt` / `fail-het-qc.txt` | FID IID of samples failing heterozygosity ±3SD filter. | Remove heterozygosity outliers before GWAS/PRS | PLINK2 (`--remove`), any text tool | `readr::read_table("het_fail_ind.txt", col_names=F)` | `pandas.read_table("het_fail_ind.txt", header=None)` |

- **R**: `readr::read_tsv("file.ext")`, `data.table::fread("file.ext")`, or `utils::read.table("file.ext", header=TRUE)`
- **Py**: `pandas.read_table("file.ext")`, `polars.read_csv("file.ext", separator="\t")`


## Completion sentinel files

| Ext | Info | Use | Tools |
|-----|------|-----|-------|
| `.done` | Empty file created on successful rule completion. | Pipeline orchestration (Snakemake `--until`), checkpoint detection | `test -f file.done` in shell; any text tool |


## See also

- [QC Pipeline Tutorial](tutorial_qc_pipeline.md) — full pipeline walkthrough
- [Ancestry Classification Tutorial](tutorial_ancestry_classification.md) — ancestry pipeline details
- [Example config](../config/example_config.yaml) — pipeline configuration reference
- [PLINK2 format specification](https://www.cog-genomics.org/plink/2.0/formats) — official PLINK2 documentation
