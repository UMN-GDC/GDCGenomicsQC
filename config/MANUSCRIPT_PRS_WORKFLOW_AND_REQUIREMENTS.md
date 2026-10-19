# Manuscript PRS Simulation Benchmark: Complete MSI README

## Purpose

This is the reproducible manuscript benchmark for six PRS methods on simulated
AFR/EUR CTSLEB genotype data:

- Single-ancestry C+T
- Single-ancestry LDpred2
- Single-ancestry lassosum2
- PRS-CSx
- CT-SLEB
- PROSPER

It creates ten bivariate phenotype replicates, then compares all methods on
the same locked AFR test participants. The primary endpoint is held-out AFR
prediction R2. The workflow also saves Slurm accounting needed to compare
runtime, core-hours, memory use, and completion/failure status.

Run the workflow only from the existing MSI project:

```text
/scratch.global/saonli/GDCGenomicsQC
```

The manuscript output root is deliberately separate from earlier exploratory
runs:

```text
/scratch.global/saonli/GDCGenomicsQC/output_ctsleb_prs_manuscript
```

Do not mix these output files with `output_ctsleb_prs` or
`output_ctsleb_prs_corrected` when preparing manuscript tables.

## Analysis Design

For each of the 10 phenotype replicates, the workflow uses:

| Component | AFR participants | EUR participants | Use |
| --- | ---: | ---: | --- |
| Discovery GWAS | 90,000 | 90,000 | Generate ancestry-specific summary statistics |
| Independent LD reference | 10,000 | 10,000 | CT-SLEB clumping/reference operations |
| Tuning set | 10,000 | Not a target evaluation set | Select all data-dependent parameters |
| Locked test set | 10,000 | Not a target evaluation set | Report primary AFR prediction results |

The simulated phenotype has requested SNP heritability `h2 = 0.40` and
cross-ancestry causal-effect correlation `rho = 0.80`. The simulation and PRS
inputs are restricted to the positional intersection of the AFR and EUR
genotypes with the PRS-CSx 1KG HapMap3 metadata. This means every method is
evaluated on the same eligible analysis panel.

The shared AFR tuning/test allocation is used by every method:

- C+T, LDpred2, and lassosum2 select their hyperparameters on tuning only.
- PRS-CSx produces a score for the locked AFR test set.
- CT-SLEB selects an empirical-Bayes score on tuning and evaluates that one
  selected score on test.
- PROSPER tunes its ensemble on tuning and evaluates it once on test.

Never select a method, threshold, or hyperparameter using test-set R2.

## What Already Exists on MSI

These large data and software resources stay on MSI. Do not upload copies of
them from the local computer.

| Requirement | MSI location |
| --- | --- |
| Project root | `/scratch.global/saonli/GDCGenomicsQC` |
| AFR genotype prefix | `/scratch.global/saonli/GDCGenomicsQC/CTSLEB/AFR/CTSLEB_AFR.{bed,bim,fam}` |
| EUR genotype prefix | `/scratch.global/saonli/GDCGenomicsQC/CTSLEB/EUR/CTSLEB_EUR.{bed,bim,fam}` |
| CT-SLEB program and reference data | `/scratch.global/saonli/GDCGenomicsQC/prs_resources/software/CTSLEB` and `/scratch.global/saonli/GDCGenomicsQC/prs_resources/ld/ct_sleb` |
| PROSPER program and downloaded LD matrices | `/scratch.global/saonli/GDCGenomicsQC/prs_resources/software/PROSPER` and `/scratch.global/saonli/GDCGenomicsQC/prs_resources/ld/prosper` |
| PRS-CSx program | `/projects/standard/gdc/public/prs_methods/scripts/PRScsx` |
| PRS-CSx 1KG reference LD | `/projects/standard/gdc/public/prs_methods/ref/ref_PRScsx/1kg_ref` |
| PLINK 1.9 | `/common/software/install/migrated/plink/1.90b6.10/bin/plink` |
| PLINK 2 | `/projects/standard/gdc/public/plink2` |
| Snakemake controller environment | `/projects/standard/gdc/public/envs/snakemake` |
| PROSPER and CT-SLEB R environment | `/scratch.global/saonli/conda-envs/prosper/bin/Rscript` |

## Files Uploaded Once From the Local Computer

The local source directory is:

```text
/Users/saonli/R25_preparing_materials/GDCGenomicsQC_fresh
```

Do not upload files one at a time. The upload helper transfers the complete,
versioned manuscript bundle in one rsync operation and creates a timestamped
backup of replaced MSI files first. The bundle includes:

- Snakemake rules: `simulatePhenotype.smk`, `preparePRSInputs.smk`,
  `prsPipelines.smk`, and `prosper_multi_prs_rule.smk`.
- Shared simulation and PRS-input scripts.
- CT-SLEB and PROSPER wrappers.
- The single-ancestry runner plus C+T, LDpred2, lassosum2, PRS-CSx, and common
  scoring/evaluation helpers.
- The manuscript configuration, controller, replicate-config generator, and
  result/compute-time collectors.
- ABCD-ready split and evaluation helpers plus the workflow documentation.

From the local computer, run exactly:

```bash
cd /Users/saonli/R25_preparing_materials/GDCGenomicsQC_fresh
bash /Users/saonli/R25_preparing_materials/GDCGenomicsQC_fresh/config/upload_ctsleb_manuscript_bundle_to_msi.sh
```

The helper uploads to `/scratch.global/saonli/GDCGenomicsQC` and saves
replaced files under:

```text
/scratch.global/saonli/GDCGenomicsQC/manuscript_prs_backups/<timestamp>/
```

The controller and configuration files that drive the benchmark are:

```text
/scratch.global/saonli/GDCGenomicsQC/config/msiCTSLEBPRS_manuscript.yaml
/scratch.global/saonli/GDCGenomicsQC/config/run_ctsleb_prs_benchmark.sh
/scratch.global/saonli/GDCGenomicsQC/config/make_ctsleb_prs_replicate_configs.sh
/scratch.global/saonli/GDCGenomicsQC/config/collect_ctsleb_prs_results.py
/scratch.global/saonli/GDCGenomicsQC/config/collect_prs_compute_metrics.py
```

## One-Time Preflight on Agate

Log in to an Agate **login node**, not a Slurm compute allocation. The
controller itself submits Slurm jobs, so do not launch it inside `srun`,
`salloc`, or an existing batch job.

```bash
ssh saonli@agate.msi.umn.edu
source /projects/standard/gdc/public/envs/load_miniconda3.sh
conda activate /projects/standard/gdc/public/envs/snakemake
cd /scratch.global/saonli/GDCGenomicsQC
```

Confirm that the required files are present before submitting anything:

```bash
ls -lh \
  /scratch.global/saonli/GDCGenomicsQC/config/msiCTSLEBPRS_manuscript.yaml \
  /scratch.global/saonli/GDCGenomicsQC/config/run_ctsleb_prs_benchmark.sh \
  /scratch.global/saonli/GDCGenomicsQC/CTSLEB/AFR/CTSLEB_AFR.bed \
  /scratch.global/saonli/GDCGenomicsQC/CTSLEB/EUR/CTSLEB_EUR.bed \
  /scratch.global/saonli/GDCGenomicsQC/prs_resources/ld/prosper/ref_bim.txt
```

Check the complete DAG without submitting jobs:

```bash
cd /scratch.global/saonli/GDCGenomicsQC/workflow
snakemake -n -p simulateBivariatePhenotypes \
  --profile /scratch.global/saonli/GDCGenomicsQC/profiles/hpc \
  --configfile /scratch.global/saonli/GDCGenomicsQC/config/msiCTSLEBPRS_manuscript.yaml
```

## Run the Complete Simulation Plan

The benchmark controller is the only command needed for the full plan. It
will generate the 10 phenotype replicates if the simulation manifest does not
already exist, generate one PRS configuration per replicate, submit the
selected methods, and finally collect prediction and compute-time results.

First run one replicate as an integration check:

```bash
cd /scratch.global/saonli/GDCGenomicsQC
N_REPLICATES=1 \
bash /scratch.global/saonli/GDCGenomicsQC/config/run_ctsleb_prs_benchmark.sh
```

This one-replicate run executes the following Snakemake targets for `rep01`:

```text
run_singleAncestryPRS
runMultiAncestryPRSCSx
runMultiAncestryPROSPER
runMultiAncestryCTSLEB
```

The single-ancestry target runs C+T, LDpred2, and lassosum2 because the
manuscript configuration specifies `single_ancestry_flags: "-c -s -l"`.

After all six methods have a complete held-out AFR metric for `rep01`, resume
the full set of ten replicates:

```bash
cd /scratch.global/saonli/GDCGenomicsQC
N_REPLICATES=10 \
bash /scratch.global/saonli/GDCGenomicsQC/config/run_ctsleb_prs_benchmark.sh
```

The controller is resumable. It reuses valid completed outputs and retries
partial jobs with `--rerun-incomplete`. Run only one copy of this controller
against `output_ctsleb_prs_manuscript` at a time.

### Optional Method-Specific Retry

Use this only after the main controller has stopped, `squeue -u saonli` is
empty, and no other Snakemake controller is running in this workflow. It is
useful for retrying a single expensive method while preserving completed ones.

```bash
cd /scratch.global/saonli/GDCGenomicsQC/workflow

snakemake -p runMultiAncestryPRSCSx \
  --profile /scratch.global/saonli/GDCGenomicsQC/profiles/hpc \
  --configfile /scratch.global/saonli/GDCGenomicsQC/config/msiCTSLEBPRS_manuscript_rep01.yaml \
  --rerun-incomplete \
  --latency-wait 180 \
  --set-resources runMultiAncestryPRSCSx:mem_mb=64000 \
  --set-resources runMultiAncestryPRSCSx:runtime=2160
```

Each `--set-resources` value must include the rule name. In particular, do
not write a bare `runtime=2160`, because Snakemake will reject it.

## Monitor Jobs and Progress

View all active jobs:

```bash
squeue -u saonli
```

The controller records every submitted Slurm job ID here:

```text
/scratch.global/saonli/GDCGenomicsQC/output_ctsleb_prs_manuscript/benchmark_compute_jobids.tsv
```

Follow each method log:

```bash
tail -f /scratch.global/saonli/GDCGenomicsQC/output_ctsleb_prs_manuscript/logs/runSingleAncestryPRS_AFR.log
tail -f /scratch.global/saonli/GDCGenomicsQC/output_ctsleb_prs_manuscript/logs/runMultiAncestryPRSCSx.log
tail -f /scratch.global/saonli/GDCGenomicsQC/output_ctsleb_prs_manuscript/logs/runMultiAncestryCTSLEB.log
tail -f /scratch.global/saonli/GDCGenomicsQC/output_ctsleb_prs_manuscript/logs/runMultiAncestryPROSPER.log
```

For a submitted Slurm job, inspect its final resource state:

```bash
sacct -j JOB_ID --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS,ReqMem
```

PRS-CSx writes paired posterior-weight files as it completes chromosomes. For
`rep01`, check their count with:

```bash
find /scratch.global/saonli/GDCGenomicsQC/output_ctsleb_prs_manuscript/prs_inputs_ctsleb/AFR_EUR/rep01/prs_pipeline/PRScsx \
  -name 'PRScsx_AFR_pst_eff_*_chr*.txt' -size +0c | wc -l
```

A value below 22 means that PRS-CSx has not completed all autosomes. A Slurm
state of `TIMEOUT` means the job ended because of the wall-time limit, not
because it necessarily ran out of memory. Retry it through the controller or
the method-specific command above with the 36-hour resource request.

## Required Outputs and Completion Checks

Set these once for convenient inspection:

```bash
OUT=/scratch.global/saonli/GDCGenomicsQC/output_ctsleb_prs_manuscript
ROOT=$OUT/prs_inputs_ctsleb/AFR_EUR/rep01
```

### Simulation Truth and Cohort Allocation

```bash
column -s $'\t' -t "$OUT/simulations/AFR_EUR/simulation_truth/simulation_manifest.tsv"
column -s $'\t' -t "$ROOT/metadata/cohort_allocation.tsv"
```

The simulation manifest must show 10 rows, realized heritability close to 0.4,
and realized cross-ancestry effect correlation close to 0.8. The allocation
file verifies the discovery GWAS, independent LD reference, tuning, and test
sample sizes.

### Method-Level Outputs

| Method | Completion marker | Held-out AFR metric |
| --- | --- | --- |
| C+T | `$ROOT/single_ancestry_AFR.done` | `$ROOT/single_ancestry_AFR/prs_pipeline/CT/C+T_heldout_metrics.tsv` |
| LDpred2 | `$ROOT/single_ancestry_AFR.done` | `$ROOT/single_ancestry_AFR/prs_pipeline/LDpred2/LDpred2-inf_heldout_metrics.tsv` and `LDpred2-grid_heldout_metrics.tsv` |
| lassosum2 | `$ROOT/single_ancestry_AFR.done` | `$ROOT/single_ancestry_AFR/prs_pipeline/lassosum2/lassosum2_heldout_metrics.tsv` |
| PRS-CSx | `$ROOT/method_runs/multi_prscsx.done` | `$ROOT/prs_pipeline/PRScsx/AFR_PRS_sscore_metrics.tsv` |
| CT-SLEB | `$ROOT/method_runs/multi_ctsleb.done` | `$ROOT/method_runs/multi_ctsleb/AFR_EUR_ctsleb_fair_metrics.tsv` |
| PROSPER | `$ROOT/method_runs/multi_prosper.done` | `$ROOT/method_runs/multi_prosper/after_ensemble_AFR/PROSPER_fair_metrics.tsv` |

For example, inspect the CT-SLEB held-out result:

```bash
column -s $'\t' -t "$ROOT/method_runs/multi_ctsleb/AFR_EUR_ctsleb_fair_metrics.tsv"
```

The selected score must be evaluated on `test_n = 10000`; its `test_r2`, not
`tuning_r2`, is the manuscript endpoint.

## Final Manuscript Tables

When all ten replicates have completed, the controller writes:

```text
/scratch.global/saonli/GDCGenomicsQC/output_ctsleb_prs_manuscript/benchmark_results_raw.tsv
/scratch.global/saonli/GDCGenomicsQC/output_ctsleb_prs_manuscript/benchmark_compute_metrics.tsv
/scratch.global/saonli/GDCGenomicsQC/output_ctsleb_prs_manuscript/benchmark_compute_jobids.tsv
/scratch.global/saonli/GDCGenomicsQC/output_ctsleb_prs_manuscript/benchmark_run_manifest.tsv
```

Review them with:

```bash
OUT=/scratch.global/saonli/GDCGenomicsQC/output_ctsleb_prs_manuscript
column -s $'\t' -t "$OUT/benchmark_results_raw.tsv" | less -S
column -s $'\t' -t "$OUT/benchmark_compute_metrics.tsv" | less -S
column -s $'\t' -t "$OUT/benchmark_run_manifest.tsv" | less -S
```

For the manuscript, summarize across the 10 paired phenotype replicates:

- Mean held-out AFR R2, standard deviation, and confidence interval per method.
- Paired differences in held-out R2 between methods using the same replicate.
- Completed, failed, timed-out, and retried run counts per method.
- Wall-clock hours, allocated core-hours, total CPU time, peak RSS, and
  requested memory from `benchmark_compute_metrics.tsv`.
- The selected method parameters and the exact score/metric source files,
  retained in the replicate directories for audit.

Do not compare methods based on tuning-set R2, package-internal validation
statistics, or incomplete chromosome-level PRS-CSx weight files.

## Safe Recovery Rules

- Do not run two controllers or manual Snakemake commands that write to the
  same output root concurrently.
- Do not use `--unlock` while a submitted workflow job or controller is still
  active. First check `squeue -u saonli` and the original terminal running the
  controller.
- If no workflow job and no controller process remains, then clear a stale lock
  from the workflow directory with:

```bash
cd /scratch.global/saonli/GDCGenomicsQC/workflow
snakemake --unlock \
  --profile /scratch.global/saonli/GDCGenomicsQC/profiles/hpc \
  --configfile /scratch.global/saonli/GDCGenomicsQC/config/msiCTSLEBPRS_manuscript_rep01.yaml
```

- Do not remove output files to force a rerun. Use the same controller command
  or `--rerun-incomplete`; valid completed outputs will be reused.

## ABCD Adaptation

The same protected discovery/tuning/testing principle applies to ABCD. Before
using this workflow structure with real ABCD data, require:

- A trait-matched external GWAS with documented non-overlap with ABCD.
- A pre-specified phenotype definition and covariate model.
- An independent, ancestry-matched PLINK LD reference.
- Harmonized target and reference variant-extract files for one locked common
  analysis panel.
- Family-level splitting, so related participants are never separated across
  tuning and test sets.
- The same tuning participants, testing participants, phenotype, covariates,
  and endpoint across all methods.

The current ABCD entry point is the single-ancestry comparison wrapper:

```text
/scratch.global/saonli/GDCGenomicsQC/workflow/scripts/run_abcd_single_ancestry_prs_pipeline.sh
```

Its regulated configuration is:

```text
/scratch.global/saonli/GDCGenomicsQC/config/abcd_single_ancestry_prs.regulated.config.sh
```

Adding PRS-CSx, CT-SLEB, or PROSPER to ABCD needs a separate multi-ancestry
run specification with ancestry-specific ABCD targets, compatible external
GWAS summary statistics, and the required method-specific reference resources.
