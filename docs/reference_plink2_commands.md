---
title: "Common PLINK2 Operations Reference"
---

## Overview

This page documents every PLINK2 command used in the GDCGenomicsQC pipeline,
organized in a single reference table. Use it to extend pipeline outputs — run
your own GRM, extract unrelated samples, compute QC metrics on custom subsets,
or apply filters with different thresholds.

```{tip}
Throughout this guide, `PREFIX` is the file path without the extension
(e.g. `/path/to/data/f1` for `f1.pgen` + `f1.pvar` + `f1.psam`).
```

## Table legend

- **Report** flags produce summary/metric files without modifying the dataset.
- **Filter** flags remove samples or variants (used with `--make-pgen` to write a
  filtered copy).
- Commands with neither label serve other purposes (I/O, LD, relatedness, etc.).

## Data format conversion

| Operation | Flag(s) | Description | Example | Output(s) |
|-----------|---------|-------------|---------|-----------|
| PLINK1 → PLINK2 | `--bfile` `--make-pgen` | Convert BED/BIM/FAM to PGEN/PVAR/PSAM | `plink2 --bfile $PREFIX --make-pgen --out ${PREFIX}_pgen` | `.pgen` `.pvar` `.psam` |
| PLINK2 → PLINK1 | `--pfile` `--make-bed` | Convert PGEN/PVAR/PSAM to BED/BIM/FAM | `plink2 --pfile $PREFIX --make-bed --out ${PREFIX}_bed` | `.bed` `.bim` `.fam` |
| PLINK2 → VCF | `--pfile` `--export vcf bgz` | Export to BGZF-compressed VCF | `plink2 --pfile $PREFIX --export vcf bgz --out $PREFIX` | `.vcf.gz` |
| VCF → PLINK2 | `--vcf` `--make-pgen` | Import VCF to PGEN. Add `--snps-only` and `--rm-dup force-first` | `plink2 --vcf input.vcf.gz --make-pgen --out $PREFIX` | `.pgen` `.pvar` `.psam` |

*Pipeline source: `convertPlink.smk` lines 98–111, 157*

## Sample and variant selection

All flags in this section are **Filter** type. Each removes a subset of the
data; combine them in a single command when multiple filters are needed.

| Operation | Flag(s) | Description | Example | Output(s) |
|-----------|---------|-------------|---------|-----------|
| Keep samples | `--keep` | Retain only samples listed in file (FID IID per line) | `plink2 --pfile $PREFIX --keep samples.txt --make-pgen --out ${PREFIX}_keep` | Updated `.pgen` `.pvar` `.psam` |
| Remove samples | `--remove` | Exclude samples listed in file (same format as --keep) | `plink2 --pfile $PREFIX --remove samples.txt --make-pgen --out ${PREFIX}_rem` | Updated `.pgen` `.pvar` `.psam` |
| Extract variants | `--extract` | Retain only variants listed in file (one ID per line) | `plink2 --pfile $PREFIX --extract vars.txt --make-pgen --out ${PREFIX}_ext` | Updated `.pgen` `.pvar` `.psam` |
| Exclude variants | `--exclude` | Remove variants listed in file | `plink2 --pfile $PREFIX --exclude vars.txt --make-pgen --out ${PREFIX}_excl` | Updated `.pgen` `.pvar` `.psam` |
| Filter by chr | `--chr` | Retain only specified chromosomes | `plink2 --pfile $PREFIX --chr 1-22 --make-pgen --out ${PREFIX}_auto` | Updated `.pgen` `.pvar` `.psam` |

```{tip}
**Combining all four:** PLINK2 applies `--remove` → `--keep` → `--exclude` →
`--extract` in that order, so you can use them together:
```
plink2 --pfile $PREFIX \
  --remove remove.txt --keep keep.txt \
  --exclude exclude.txt --extract extract.txt \
  --make-pgen --out ${PREFIX}_filtered
```
```

*Pipeline source: `convertPlink.smk` lines 116–139, `filterStandard.sh` lines 25–38, `Standard_QC.smk` lines 46–48*

## Missingness

| Operation | Flag(s) | Type | Description | Example | Output(s) |
|-----------|---------|------|-------------|---------|-----------|
| Compute missingness | `--missing` | Report | Per-sample and per-variant missing call rates | `plink2 --pfile $PREFIX --missing --out ${PREFIX}_qc` | `.smiss` `.vmiss` |
| Filter variants by missingness | `--geno` | Filter | Remove variants with call rate below threshold (e.g. 0.02 = 2%) | `plink2 --pfile $PREFIX --geno 0.02 --make-pgen --out ${PREFIX}_geno002` | Updated `.pgen` `.pvar` `.psam` |
| Filter samples by missingness | `--mind` | Filter | Remove samples with genotype rate below threshold (e.g. 0.02 = 2%) | `plink2 --pfile $PREFIX --mind 0.02 --make-pgen --out ${PREFIX}_mind002` | Updated `.pgen` `.pvar` `.psam` |

*Pipeline source: `initialFilter.sh` lines 18–24, `convertPlink.smk` line 155*

## Allele frequency

| Operation | Flag(s) | Type | Description | Example | Output(s) |
|-----------|---------|------|-------------|---------|-----------|
| Compute frequencies | `--freq` | Report | Allele frequencies for all variants | `plink2 --pfile $PREFIX --freq --out $PREFIX` | `.afreq` |
| Frequency by ancestry | `--keep` `--freq` | Report | Frequencies restricted to a sample subset | `plink2 --pfile $PREFIX --keep keep_EUR.txt --freq --out ${PREFIX}_EUR_freq` | `.afreq` |
| Filter by MAF | `--maf` | Filter | Remove variants with MAF below threshold (e.g. 0.01) | `plink2 --pfile $PREFIX --maf 0.01 --make-pgen --out ${PREFIX}_maf01` | Updated `.pgen` `.pvar` `.psam` |

*Pipeline source: `initialFilter.sh` line 30, `filterStandard.sh` lines 13–14, `Standard_QC.smk` line 53*

## Hardy-Weinberg equilibrium

| Operation | Flag(s) | Type | Description | Example | Output(s) |
|-----------|---------|------|-------------|---------|-----------|
| Compute HWE p-values | `--hardy` | Report | HWE exact test per variant (founders only by default) | `plink2 --pfile $PREFIX --hardy --out $PREFIX` | `.hardy` |
| Filter by HWE | `--hwe` | Filter | Remove variants with HWE p below threshold | `plink2 --pfile $PREFIX --hwe 1e-6 --make-pgen --out ${PREFIX}_hwe6` | Updated `.pgen` `.pvar` `.psam` |
| Filter by HWE (scaled) | `--hwe p k` | Filter | Sample-size-scaled threshold: `p × 10^(-n×k)`. Greer et al. (2024) recommends `k=0.001` | `plink2 --pfile $PREFIX --hwe 1e-6 0.001 --make-pgen --out ${PREFIX}_hwe_scaled` | Updated `.pgen` `.pvar` `.psam` |
| Extract low-p HWE (for review) | `awk` on `.hardy` | Diagnostic | Write variants with p < 1e-5 to a separate file (not a filter) | `awk '$9 < 1e-5' $PREFIX.hardy > zoomhwe.hwe` | `.hwe` (custom) |

*Pipeline source: `initialFilter.sh` lines 31, `filterStandard.sh` lines 17–20, `Standard_QC.smk` lines 55–61*

## Heterozygosity / inbreeding

| Operation | Flag(s) | Type | Description | Example | Output(s) |
|-----------|---------|------|-------------|---------|-----------|
| LD pruning for het | `--indep-pairwise` | LD | Generate LD-pruned variant list (50/5/0.2) | `plink2 --pfile $PREFIX --indep-pairwise 50 5 0.2 --out ${PREFIX}_het` | `.prune.in` `.prune.out` |
| Compute heterozygosity | `--extract` `--het` | Report | Per-sample inbreeding coefficient F on LD-pruned variants | `plink2 --pfile $PREFIX --extract ${PREFIX}_het.prune.in --het --out ${PREFIX}_het` | `.het` |
| Remove het outliers | `--remove` | Filter | Remove samples with F beyond ±3SD from mean | `plink2 --pfile $PREFIX --remove het_fail_ind.txt --make-pgen --out ${PREFIX}_nohetout` | Updated `.pgen` `.pvar` `.psam` |

The outlier detection requires an R step first:
```r
het <- read.table("${PREFIX}_het.het", header=TRUE)
het_fail <- subset(het, F < mean(F) - 3*sd(F) | F > mean(F) + 3*sd(F))
write.table(het_fail[, c("FID", "IID")], "het_fail_ind.txt", row.names=FALSE)
```

*Pipeline source: `initialFilter.sh` lines 32–33, `filterStandard.sh` lines 27–40, `heterozygosity_outliers_list.R`*

## Sex check

| Operation | Flag(s) | Type | Description | Example | Output(s) |
|-----------|---------|------|-------------|---------|-----------|
| Check sex | `--check-sex` | Report | Infer sex from X-chromosome het; compare to reported sex | `plink2 --pfile $PREFIX --check-sex --out $PREFIX` | `.sexcheck` |
| Remove discrepant | `--remove` | Filter | Remove samples where STATUS == PROBLEM | `grep PROBLEM $PREFIX.sexcheck | awk '{print \$1,\$2}' > sex_disc.txt; plink2 --pfile $PREFIX --remove sex_disc.txt --make-pgen --out ${PREFIX}_sexok` | Updated `.pgen` `.pvar` `.psam` |

*Pipeline source: `Standard_QC.smk` lines 46–48, 128–130*

## LD pruning

| Operation | Flag(s) | Type | Description | Example | Output(s) |
|-----------|---------|------|-------------|---------|-----------|
| LD pruning (standard) | `--indep-pairwise` | LD | Window 500 kb, step 10, r² 0.1. For PCA, heritability, GRM | `plink2 --pfile $PREFIX --indep-pairwise 500 10 0.1 --out $PREFIX` | `.prune.in` `.prune.out` |
| Apply prune list | `--extract` | Filter | Keep only LD-pruned variants | `plink2 --pfile $PREFIX --extract $PREFIX.prune.in --make-pgen --out ${PREFIX}_ldp` | Updated `.pgen` `.pvar` `.psam` |
| LD pruning (het) | `--indep-pairwise` | LD | Window 50 vars, step 5, r² 0.2. For heterozygosity | `plink2 --pfile $PREFIX --indep-pairwise 50 5 0.2 --out ${PREFIX}_het` | `.prune.in` `.prune.out` |

*Pipeline source: `initialFilter.sh` lines 26–27, 32, `filterStandard.sh` lines 27, 44–45, `Standard_QC.smk` line 63*

## Relatedness

| Operation | Flag(s) | Type | Description | Example | Output(s) |
|-----------|---------|------|-------------|---------|-----------|
| KING kinship | `--make-king` | Relatedness | Pairwise KING-robust kinship (LD-aware, no pruning needed) | `plink2 --pfile $PREFIX --make-king --out $PREFIX` | `.kin0` |
| KING triangle | `--make-king triangle` | Relatedness | Lower-triangular KING output (like GCTA GRM layout) | `plink2 --pfile $PREFIX --make-king triangle --out $PREFIX` | `.kin` |
| Extract unrelated | `--king-cutoff` | Filter | Remove one per pair with kinship > 0.0884 (3rd-degree) | `plink2 --pfile $PREFIX --king-cutoff 0.0884 --make-pgen --out ${PREFIX}_unrel` | `.king.cutoff.in.id` `.king.cutoff.out.id` + updated `.pgen` |
| Extract unrelated (manual) | `--keep` | Filter | Keep only samples in pre-computed unrelated list | `plink2 --pfile $PREFIX --keep unrelated.id --make-pgen --out ${PREFIX}_unrel` | Updated `.pgen` `.pvar` `.psam` |
| GRM | `--make-grm-bin` | Relatedness | Lower-triangular genetic relationship matrix (float32) | `plink2 --pfile $PREFIX --make-grm-bin --out $PREFIX` | `.grm.bin` `.grm.id` `.grm.N.bin` |

```{tip}
**GRM for heritability — full workflow:**
```bash
plink2 --pfile $PREFIX --indep-pairwise 500 10 0.1 --out $PREFIX        # LD prune
plink2 --pfile $PREFIX --keep unrelated.txt --make-pgen --out ${PREFIX}_u   # unrelated only
plink2 --pfile ${PREFIX}_u --extract $PREFIX.prune.in --make-pgen --out ${PREFIX}_u_ldp  # apply LD
plink2 --pfile ${PREFIX}_u_ldp --make-grm-bin --out ${PREFIX}_grm           # build GRM
```
```

*Pipeline source: `Relatedness.smk` lines 41–69, `simulatePhenotype2.smk` line 109*

## PCA

| Operation | Flag(s) | Type | Description | Example | Output(s) |
|-----------|---------|------|-------------|---------|-----------|
| Compute PCs | `--pca N` | Dim. reduction | N principal components + eigenvalues | `plink2 --pfile $PREFIX --pca 20 --out $PREFIX` | `.eigenvec` `.eigenval` |
| PCA with loadings | `--pca N allele-wts` | Dim. reduction | Also output per-SNP loading weights | `plink2 --pfile $PREFIX --pca 20 allele-wts --out $PREFIX` | + `.eigenvec.allele` |
| Project new samples | `--read-freq` `--score` | Dim. reduction | Score study samples using existing PC loadings from reference | `plink2 --pfile $STUDY --read-freq $REF.acount --score $REF.eigenvec.allele header iid-read 2 3 4 --out ${STUDY}_proj` | `.sscore` |

*Pipeline source: `PCAreference.smk` lines 107, 119, `simulatePhenotype2.smk` line 133*

## Allele alignment

| Operation | Flag(s) | Type | Description | Example | Output(s) |
|-----------|---------|------|-------------|---------|-----------|
| Align to reference | `--fa` `--ref-from-fa force` | Alignment | Force REF/ALT to match reference genome FASTA (auto strand-flip) | `plink2 --pfile $PREFIX --fa ref.fa --ref-from-fa force --make-pgen --out ${PREFIX}_aligned` | Updated `.pgen` `.pvar` `.psam` |
| Rename variant IDs | `--set-all-var-ids` | Alignment | Set IDs to CHR:POS:REF:ALT format | `plink2 --pfile $PREFIX --set-all-var-ids 'chr@:#:\$r:\$a' --make-pgen --out ${PREFIX}_id` | Updated `.pvar` |
| Strand flip | `--flip` | Alignment | Flip strand of specified variants | `plink2 --pfile $PREFIX --flip flip_list.txt --make-pgen --out ${PREFIX}_flipped` | Updated `.pgen` `.pvar` `.psam` |

*Pipeline source: `convertPlink.smk` lines 157–204, 343–373*

## Merging

| Operation | Flag(s) | Type | Description | Example | Output(s) |
|-----------|---------|------|-------------|---------|-----------|
| Concatenate per-chr | `--pmerge-list` | Merge | Merge per-chromosome PGENs into single genome-wide set | `echo /path/to/chr1/f1 > mergelist.txt; echo /path/to/chr2/f1 >> mergelist.txt; plink2 --pmerge-list mergelist.txt --make-pgen --out ${PREFIX}_merged` | `.pgen` `.pvar` `.psam` |

*Pipeline source: `convertPlink.smk` line 401*

## Complete QC pipeline (command-line equivalent)

```bash
INPUT="f1"

# --- Initial QC ---
plink2 --pfile $INPUT --mind 0.1 --make-pgen --out ${INPUT}_2
plink2 --pfile ${INPUT}_2 --geno 0.02 --make-pgen --out ${INPUT}_3
plink2 --pfile ${INPUT}_3 --mind 0.02 --make-pgen --out ${INPUT}_qc

# Metrics
plink2 --pfile ${INPUT}_qc --freq --out ${INPUT}_qc
plink2 --pfile ${INPUT}_qc --hardy --out ${INPUT}_qc
plink2 --pfile ${INPUT}_qc --indep-pairwise 50 5 0.2 --out ${INPUT}_het
plink2 --pfile ${INPUT}_qc --extract ${INPUT}_het.prune.in --het --out ${INPUT}_het

# --- Standard QC ---
plink2 --pfile ${INPUT}_qc --maf 0.01 --make-pgen --out ${INPUT}_s2
plink2 --pfile ${INPUT}_s2 --hardy --out ${INPUT}_s2
plink2 --pfile ${INPUT}_s2 --hwe 1e-6 --make-pgen --out ${INPUT}_s3a
plink2 --pfile ${INPUT}_s3a --hwe 1e-10 --make-pgen --out ${INPUT}_s3

# Heterozygosity outliers
plink2 --pfile ${INPUT}_s3 --indep-pairwise 50 5 0.2 --out ${INPUT}_indep
plink2 --pfile ${INPUT}_s3 --extract ${INPUT}_indep.prune.in --het --out ${INPUT}_check
# R: het <- read.table("${INPUT}_check.het", header=TRUE)
#     het_fail <- subset(het, F < mean(F)-3*sd(F) | F > mean(F)+3*sd(F))
#     write.table(het_fail[,c("FID","IID")], "het_fail.txt", row.names=FALSE)
plink2 --pfile ${INPUT}_s3 --remove het_fail.txt --make-pgen --out ${INPUT}_s4

# Final LD-pruned output
plink2 --pfile ${INPUT}_s4 --indep-pairwise 500 10 0.1 --out ${INPUT}_final
plink2 --pfile ${INPUT}_s4 --extract ${INPUT}_final.prune.in --make-pgen --out ${INPUT}_final.LDpruned
```

## PLINK2 flag quick-reference

| Flag | Category | Purpose |
|------|----------|---------|
| `--pfile` / `--bfile` | I/O | Open PLINK2 / PLINK1 dataset |
| `--make-pgen` / `--make-bed` | I/O | Write PLINK2 / PLINK1 dataset |
| `--out` | I/O | Output prefix |
| `--keep` / `--remove` | Filter | Sample inclusion / exclusion |
| `--extract` / `--exclude` | Filter | Variant inclusion / exclusion |
| `--chr` | Filter | Chromosome range |
| `--maf` | Filter | Minor allele frequency threshold |
| `--geno` | Filter | Variant missingness threshold |
| `--mind` | Filter | Sample missingness threshold |
| `--hwe` | Filter | HWE p-value threshold |
| `--snps-only` | Filter | Restrict to SNPs |
| `--rm-dup` | Filter | Remove duplicate variant IDs |
| `--missing` | Report | Sample and variant missingness |
| `--freq` | Report | Allele frequencies |
| `--hardy` | Report | HWE exact test p-values |
| `--het` | Report | Heterozygosity / inbreeding F |
| `--check-sex` | Report | Sex inference from X-chromosome |
| `--indep-pairwise` | LD | LD pruning |
| `--make-king` | Relatedness | KING-robust kinship |
| `--king-cutoff` | Relatedness | Unrelated sample extraction |
| `--make-grm-bin` | Relatedness | Genetic relationship matrix |
| `--pca` | Dim. reduction | Principal components |
| `--score` | Dim. reduction | Polygenic scoring / projection |
| `--flip` | Alignment | Strand flip |
| `--ref-from-fa` | Alignment | REF allele from FASTA |
| `--fa` | Alignment | Reference genome FASTA |
| `--set-all-var-ids` | Alignment | Rename variant IDs |
| `--pmerge-list` | Merge | Concatenate per-chr PGENs |
| `--threads` | Performance | CPU thread count |
| `--memory` | Performance | RAM limit (MB) |

## See also

- [Pipeline Output File Reference](reference_output_formats.md) — detailed
  descriptions of every output format
- [QC Pipeline Tutorial](tutorial_qc_pipeline.md) — step-by-step pipeline
  walkthrough
- [PLINK2 documentation](https://www.cog-genomics.org/plink/2.0/) — official
  reference for all flags
