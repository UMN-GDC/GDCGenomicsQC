# CTSLEB Manuscript PRS Benchmark

## Purpose

This manuscript workflow runs the corrected AFR/EUR phenotype simulation and a
fair PRS comparison for C+T, LDpred2, lassosum2, PRS-CSx, CT-SLEB, and
PROSPER. It uses a single protected AFR tuning/test split per simulation
replicate. No test outcomes are used for hyperparameter selection.

The existing MSI project directory is:

`/scratch.global/saonli/GDCGenomicsQC`

All outputs are isolated below:

`/scratch.global/saonli/GDCGenomicsQC/output_ctsleb_prs_manuscript`

The input CTSLEB genotypes, shared PRS software, and reference panels remain
in their existing locations. The upload script backs up each replaced workflow
file under `manuscript_prs_backups/` before updating it.

## Analysis Design

- Ten bivariate phenotype replicates: requested h2=0.40 and cross-ancestry effect correlation rho=0.80.
- Each ancestry contributes a 90,000-participant discovery GWAS.
- Each ancestry has a distinct 10,000-person LD-reference panel. CT-SLEB uses
  those panels for its target and auxiliary LD-dependent clumping steps.
- The remaining AFR participants are split into a 10,000-person tuning cohort
  and a locked 10,000-person test cohort.
- Phenotypes and all analysis genotype subsets are restricted to the positional
  intersection of AFR, EUR, and the PRS-CSx 1KG HapMap3 reference metadata.
- C+T, LDpred2, and lassosum2 select parameters on tuning only.
- PRS-CSx scores the locked AFR test cohort only.
- CT-SLEB and PROSPER use method-specific files derived from that same split.
- The primary endpoint is held-out AFR R2. The truth manifest additionally verifies realized h2 and cross-ancestry effect correlation.

## Run on Agate

From the local `GDCGenomicsQC_fresh` folder, update the manuscript files once:

```bash
bash config/upload_ctsleb_manuscript_bundle_to_msi.sh
```

Run the controller only from an Agate login node:

```bash
ssh saonli@agate.msi.umn.edu
source /projects/standard/gdc/public/envs/load_miniconda3.sh
conda activate /projects/standard/gdc/public/envs/snakemake
cd /scratch.global/saonli/GDCGenomicsQC
```

First run one complete replicate as an integration test:

```bash
N_REPLICATES=1 bash config/run_ctsleb_prs_benchmark.sh
```

After every method completes and its held-out metrics pass review, resume the
full benchmark. Completed outputs are reused automatically:

```bash
N_REPLICATES=10 bash config/run_ctsleb_prs_benchmark.sh
```

To run a narrower diagnostic without one expensive method:

```bash
RUN_CTSLEB=0 N_REPLICATES=1 bash config/run_ctsleb_prs_benchmark.sh
```

Do not run two copies of the driver against this output directory at once.

## Required Outputs

Prediction results:

`output_ctsleb_prs_manuscript/benchmark_results_raw.tsv`

Slurm resource accounting:

`output_ctsleb_prs_manuscript/benchmark_compute_metrics.tsv`

Submitted-job inventory:

`output_ctsleb_prs_manuscript/benchmark_compute_jobids.tsv`

Controller logs:

`output_ctsleb_prs_manuscript/benchmark_controller_logs/`

The compute table stores Slurm job state, exit code, elapsed wall time,
allocated CPUs, allocated core-hours, total CPU time, maximum RSS, requested
memory, start/end times, and the associated Snakemake rule. These fields allow
method-level runtime and memory comparisons, including failed runs.

## Monitoring

```bash
squeue -u saonli
tail -n 100 output_ctsleb_prs_manuscript/benchmark_run_manifest.tsv
tail -n 100 output_ctsleb_prs_manuscript/logs/runMultiAncestryPRSCSx.log
```

After a run completes:

```bash
column -s $'\t' -t output_ctsleb_prs_manuscript/benchmark_results_raw.tsv | less -S
column -s $'\t' -t output_ctsleb_prs_manuscript/benchmark_compute_metrics.tsv | less -S
```

## ABCD Adaptation

The same protected-split and common-evaluation design is documented in
`config/README_prs_manuscript_workflow.md`. The ABCD runner creates a split
manifest, retains a test-score file for audit, and supports covariate-adjusted
continuous or binary endpoints. Use an external, non-overlapping discovery
GWAS and keep families together whenever related participants are retained.
