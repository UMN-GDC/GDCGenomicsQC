# Reproducibility utilities

This directory contains public synthetic fixtures, expected-output checks, and
MSI validation utilities. It contains no participant-level data.

## Portable smoke test

From a clean clone with Snakemake, conda/mamba, and the required container
runtime available:

```bash
./reproducibility/scripts/run_smoke_test.sh
```

The command builds a deterministic 60-sample, 6-variant VCF fixture; renders a
configuration with absolute paths; runs `run_initialQC`; and verifies the
required output files and minimum retained sample/variant counts. It exercises
the v1 single-file input, GRCh38 ID normalization, basic QC, and standard-QC
paths. It does not test liftover, reference projection, VAE, RFMix, Gnomix,
GWAS, release filtering, or heritability. The expected contract is in
`expected_outputs.json`, and ground truth for the fixture is in
`fixtures/truth.json`.

## MSI 1kgSynthetic global-ancestry integration validation

The MSI shared practice dataset `1kgSynthetic` contains 150,000 synthetic
samples (30,000 each AFR, AMR, EAS, EUR, and SAS) and 97,915 variants. The job
below exercises the v1 PLINK BED input, QC, reference-projected
PCA, and PCA-based ancestry-classification path with a fixed seed. It writes a
confusion matrix and accuracy/uncertainty summary against the dataset's labels.

```bash
sbatch reproducibility/jobs/msi_1kgsynthetic_global_ancestry.SLURM
```

This is a labelled integration test. The job preserves its rendered config,
evaluation summary, confusion matrix, master log, and output checksums under
the selected scratch directory.
