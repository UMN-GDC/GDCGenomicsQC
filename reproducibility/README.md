# Reproducibility package

This directory is the submission-facing contract for the software paper. It
contains no participant-level data. The small fixture is synthetic and exists
only to test workflow behavior, not to evaluate population-genetic methods.

## What reviewers can run

From a clean clone with Snakemake, conda/mamba, and the required container
runtime available:

```bash
./reproducibility/scripts/run_smoke_test.sh
```

The command builds a deterministic 12-sample, 6-variant VCF fixture; renders a
configuration with absolute paths; runs `run_initialQC`; and verifies the
required output files and minimum retained sample/variant counts. The expected
contract is in `expected_outputs.json`, and ground truth for the fixture is in
`fixtures/truth.json`.

## Before manuscript submission

1. Run the smoke test on the tagged release using both the local-conda and the
   intended sandbox/container profile. Save the logs and output manifest.
2. Add an integration fixture based on a public 1000 Genomes-derived subset for
   reference projection and local ancestry. Keep its licensing/provenance file
   with the fixture; do not present the tiny smoke test as ancestry validation.
3. Add simulated-truth benchmarks for classification, local ancestry,
   heritability, and each advertised PRS module. Each benchmark needs its input,
   seed, config, expected summary, tolerance, and evaluator.
4. Tag the exact software release and run
   `./reproducibility/scripts/capture_release_lockfiles.sh`. Archive the
   resulting conda explicit locks, container image checksums, source commit, and
   source-status files with the GigaDB/Zenodo deposit.
5. Replace all mutable `latest` container tags in the release branch with image
   digests or archive the pulled SIF images and their checksums. A tag such as
   `latest` is not a reproducible release identifier.

## Release manifest

The final archive must include: source commit and tag; clean working-tree
status; fixture checksums; all rendered configs; test commands and logs;
container digests/SIF checksums; conda explicit lockfiles; and the expected
outputs used to declare a run successful.

Do not claim end-to-end validation of a pipeline branch until its fixture,
expected result, and evaluator are deposited here.
