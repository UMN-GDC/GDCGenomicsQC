# Manuscript PRS Workflow: Simulations and ABCD

This workflow uses one design for simulated and real-data PRS benchmarks:

1. **Discovery:** PRS weights are estimated from external GWAS summary statistics. ABCD participants must not overlap a discovery GWAS.
2. **Tuning:** One protected target-cohort partition selects every data-dependent setting, including C+T thresholds, LDpred2/lassosum2 grids, CT-SLEB candidates, and PROSPER ensembles.
3. **Testing:** A disjoint locked partition is scored once per selected model. Only this partition contributes the primary manuscript endpoint.

## Reusable Cohort Contract

Every analysis must provide headered, tab-separated files with `FID`, `IID`, and `PHENO`:

- `<label>_tuning.pheno.tsv` and `<label>_testing.pheno.tsv`
- `<label>_tuning.keep.txt` and `<label>_testing.keep.txt`
- `<label>_split_assignments.tsv` and `<label>_split_manifest.json`

`workflow/scripts/prepare_prs_tuning_test_split.py` creates these files deterministically. For cohorts containing relatives, supply a headered `FID`, `IID`, `GROUP` table through `--group-file`; all participants in a family/group are assigned together. A prior unrelated-only target may omit it.

## CTSLEB Simulation Benchmark

The corrected CTSLEB configuration already produces one shared AFR tuning/test split for each replicate. All methods must consume that shared split:

- C+T, LDpred2, and lassosum2 tune on the tuning partition and write test PRS scores.
- PRS-CSx scores only the held-out target participants.
- CT-SLEB and PROSPER use the same tuning/test participants through their method-specific input formats.

The simulation truth manifest verifies the data-generating model. Primary performance is held-out prediction R2. Across replicates, report mean, standard deviation, confidence intervals, paired method differences, runtime, and convergence/failure counts.

## ABCD Analysis

Use `config/abcd_single_ancestry_prs.regulated.config.sh` with the ABCD entry point. It now creates a protected split at:

`OUTPUT_ROOT/<ancestry>/<trait>/inputs/split/`

Submit the ABCD wrapper itself through Slurm with `RUN_WITH_SBATCH=false`; this
prevents a Slurm job from submitting a second nested Slurm job. Set
`RUN_WITH_SBATCH=true` only when launching the wrapper directly from a login
node.

Before a final run:

- Replace the placeholder GWAS with a trait-matched external GWAS and document its sample size, ancestry composition, genome build, and non-overlap with ABCD.
- Supply an independent ancestry-matched PLINK LD reference through `LD_REFERENCE_BFILE_PREFIX`. The final C+T analysis uses this panel for clumping, rather than an additional split of ABCD.
- Create `TARGET_VARIANT_EXTRACT` and `LD_REFERENCE_VARIANT_EXTRACT` after harmonizing the ABCD target and external reference to the same build, alleles, and locked common variant panel. Each file is one BIM ID per line; they may use different IDs despite representing the same loci.
- For ABCD analyses retaining relatives, provide `FAMILY_GROUP_FILE`; do not split a family across tuning and testing.
- Populate `COVARIATE_FILE` and `COVARIATE_COLUMNS` with pre-specified variables, for example age, sex, recruitment site, genotyping batch, and ancestry PCs. Use only variables measured independently of the outcome.
- Run every PRS method with the same target ancestry, tuning participants, testing participants, phenotype definition, and covariate model.

Set `MANUSCRIPT_MODE=true` only after replacing the demonstration GWAS,
setting `GWAS_TRAIT_LABEL`, and documenting non-overlap through
`GWAS_TARGET_OVERLAP_CONFIRMED=true`. The generated `analysis_manifest.tsv`
records these choices, sample counts, variant counts, reference-panel details,
and the locked split beside the method outputs. Manuscript mode refuses to run
without an external LD reference and both harmonized variant-extract files.

`workflow/scripts/prs_pipeline/src/evaluate_heldout_prs.R` is the common endpoint calculator. It consumes an individual-level held-out PRS score, held-out phenotype, and optional covariates. For continuous outcomes, its primary endpoint is covariate-adjusted incremental R2. For binary outcomes, it reports incremental deviance R2 and AUC. It also writes the exact analyzed scores for audit and calibration plots.

The ABCD single-ancestry runner calls this evaluator automatically for C+T,
LDpred2, and lassosum2 after their tuning-selected test scores are produced.
C+T uses the configured independent LD reference. The current LDpred2 and
lassosum2 wrappers estimate LD from the ancestry-filtered ABCD target sample;
this should be reported as an implementation limitation unless those wrappers
are separately extended to score from an external reference.

## Required Reporting

For every method, ancestry, trait, and split report:

- discovery GWAS source and overlap assessment;
- tuning and test participant/group counts;
- selected hyperparameters, chosen using tuning data only;
- held-out sample size, effect per PRS SD, confidence interval, P value, and the common primary metric;
- score coverage, excluded participants, runtime, and any failed run;
- an ancestry-stratified comparison table. Do not pool ancestries as a substitute for stratified reporting.

Keep the split manifest and method-level score files with the manuscript outputs. They make the final comparisons reproducible and protect against accidental retuning on the held-out cohort.
