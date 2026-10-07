# Tutorial: Release Filter Pipeline

This tutorial walks through running the **Release Filter** pipeline to produce a de-identified, filtered release dataset from post-QC genomic data.

## Prerequisites

- Completed QC pipeline through at least the relatedness/ancestry stages
- A crosswalk file mapping raw IDs (pscid) to release IDs (release_candid + C/M suffix)
- Exclusion lists (optional but recommended)
- PLINK2 bed/bim/fam or pgen/pvar/psam from the QC pipeline
- Any derivative files you want filtered (GRM, eigenvecs, CNV, etc.)

## Step 1: Prepare Configuration

Create a config file (e.g., `config_release.yaml`) with the releaseFilter block:

```yaml
releaseFilter:
    enabled: true

    # Core inputs (required)
    source_fam: "/path/to/QC_passed.fam"           # All post-QC subjects (raw IDs)
    identifiers: "/path/to/release_identifiers.csv" # pscid -> release_candid (+ suffix)
    keep_list: "/path/to/keep_list.txt"            # Output for plink2 --keep
    temp_fam: "/path/to/temp.fam"                  # De-identified .fam (all QC-passing)
    crosswalk: "/path/to/release_identifiers.csv"  # For de-identification step

    # Optional exclusion sources
    exclusion_sources:
        - "/path/to/HBCDexclusions.csv"
    par_visit: "/path/to/par_visit.csv"            # Eligible subjects (raw IDs)
    batch_info: "/path/to/batch_info.txt"          # Optional output
    removed_individuals: "/path/to/removed.txt"    # Optional output

    # ID column names
    id_col: "IID"
    release_id_col: "release_candid"
    suffix_col: "relationship"                     # C/M column
    pscid_col: "pscid"

    # Derivative files to process
    derivatives:
        bed: "/path/to/data.bed"                   # Filter with plink2 --keep
        bim: "/path/to/data.bim"
        fam: "/path/to/data.fam"
        # OR pgen/pvar/psam
        # pgen: "/path/to/data.pgen"
        # pvar: "/path/to/data.pvar"
        # psam: "/path/to/data.psam"
        grm_bin: "/path/to/grm_prefix"             # .grm.bin/.grm.id/.grm.N.bin
        grm_gz: "/path/to/grm.gz"                  # Filter .grm.gz
        eigenvectors:
            - "/path/to/internal_pca_plink2.eigenvec"
        cnv_files:
            - "/path/to/CNV_slim.txt"
        cnv_col: "sample_id"
        deidentify_files:                          # Raw ID -> de-ID
            - "/path/to/external_CNV.txt"

    validate: true
```

## Step 2: Run the Pipeline

### Option A: Using the gdcgenomicsqc wrapper (MSI/Sandbox)

```bash
# Load environment
module use /projects/standard/gdc/public/GDCGenomicsQC/envs
module load gdcgenomicsMSI
conda activate snakemake

# Run
gdcgenomicsqc run_releaseFilter --configfile config_release.yaml
```

### Option B: Using snakemake directly

```bash
conda activate snakemake
cd GDCGenomicsQC/workflow
snakemake --profile=../profiles/hpc --configfile ../config_release.yaml run_releaseFilter
```

### Option C: Local testing (interactive profile)

```bash
snakemake --profile=../profiles/interactive --configfile ../config_release.yaml run_releaseFilter
```

## Step 3: Monitor Progress

The pipeline runs these stages in order:

| Stage | Rule | Description |
|-------|------|-------------|
| 1 | `buildReleaseKeep` | Build keep-list + temp.fam from QC-passed subjects |
| 2 | `filterGenomicsWithPlink2` / `filterPgenWithPlink2` | PLINK2 `--keep` on bed/bim/fam or pgen/pvar/psam |
| 3 | `deidentifyFiles` | De-identify raw ID files using crosswalk |
| 4 | `filterDerivatives` | Subset GRM, eigenvecs, CNV, generic files to keep-list |
| 5 | `validateRelease` | Verify no excluded IDs leaked, all IIDs match pattern |

Watch the SLURM logs:
```bash
# Tail the main log
tail -f logs/run_releaseFilter.log

# Or check individual rule logs
tail -f logs/buildReleaseKeep.log
```

## Step 4: Verify Outputs

After completion, check the output files:

```bash
# Keep-list
head keep_list.txt
# Expected: one release IID per line (e.g., 1234567890C)

# Temp FAM
head temp.fam
# Expected: FID=release_candid, IID=release_candid+C/M

# Filtered PLINK
ls -la *_filtered.bed *_filtered.bim *_filtered.fam

# Filtered GRM
ls -la *_filtered.grm.bin *_filtered.grm.id *_filtered.grm.N.bin

# Filtered eigenvecs
head internal_pca_plink2_filtered.eigenvec

# Validation log (should show all checks passing)
cat logs/validateRelease.log
```

## Step 5: Validation

The `validateRelease` rule runs automatically at the end and checks:

1. **Pattern check**: All IIDs match `^\d{10}[CM]$`
2. **Exclusion leak check**: No excluded IIDs appear in outputs
3. **Keep-list subset check**: All output IIDs are in the keep-list
4. **GRM dimension check**: GRM binary dimensions match temp.fam row count

If validation fails, the pipeline exits with an error and details in `logs/validateRelease.log`.

## Common Issues

### Missing input files
```
MissingInputException: Missing input files for rule buildReleaseKeep:
    /path/to/release_identifiers.csv
    /path/to/QC_passed.fam
```
**Fix**: Ensure all input paths in config are correct and files exist.

### No subjects in keep-list
```
Final release subjects (raw): 0
```
**Fix**: Check that exclusion lists and par-visit filters aren't removing all subjects. Verify crosswalk has mappings for your raw IDs.

### Validation fails: excluded IIDs leaked
```
[FAIL] keep_list.txt: 3 excluded IIDs found!
```
**Fix**: Check that exclusion sources are correctly formatted and crosswalk mappings are correct. Re-run with `snakemake --forcerun validateRelease --configfile config_release.yaml`.

### PLINK2 --keep fails
```
Error: Invalid sample ID in --keep file
```
**Fix**: Ensure keep_list.txt has one IID per line, no headers, no extra whitespace.

## Applying at Different Pipeline Stages

The release filter is **generic** and can be applied at any stage by pointing `source_fam` and `derivatives` at the appropriate intermediate outputs:

| Stage | source_fam | derivatives.bed/fam | derivatives.grm_bin | derivatives.eigenvectors |
|-------|------------|---------------------|---------------------|--------------------------|
| Post-QC | `full/f1.b38.f2.fam` | `full/f1.b38.f2.bed` | `full/f1.b38.ldpruned` | `full/internal_pca_plink2.eigenvec` |
| Post-PCA | `EUR/f1.b38.f2.fam` | `EUR/f1.b38.f2.bed` | `EUR/f1.b38.ldpruned` | `EUR/internal_pca_plink2.eigenvec` |
| Post-imputation | `full/f1.b38.f2.fam` | `full/f1.b38.f2.pgen` | `full/f1.b38.ldpruned` | `full/internal_pca_plink2.eigenvec` |

Simply update the config paths to point to the desired stage's outputs.

## Output Structure

After successful run, the release directory will contain:

```
release/
├── keep_list.txt              # plink2 --keep input
├── temp.fam                   # de-identified FAM (all QC-passing)
├── batch_info.txt             # optional batch info
├── removed_individuals.txt    # excluded IIDs
├── data_filtered.bed/.bim/.fam   # filtered PLINK
├── grm_prefix_filtered.grm.bin/.id/.N.bin  # filtered GRM
├── internal_pca_plink2_filtered.eigenvec   # filtered PCs
├── CNV_slim_filtered.txt      # filtered CNV
└── external_CNV_deid.txt      # de-identified CNV
```

## Next Steps

- Use `keep_list.txt` with `plink2 --keep` for any downstream analysis
- Use `temp.fam` with `plink2 --bfile ... --fam temp.fam` for correct .bed reading
- Share de-identified derivatives (`*_filtered.*`, `*_deid`) with collaborators
- Run GWAS/PRS/heritability on the filtered release dataset