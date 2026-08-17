# building containers
Use the .SLURM scripts

# Upload to GHCR
Register then push

```
apptainer registry login --username coffm049 docker://ghcr.io


# redirect so it doesn't take up all space
export APPTAINER_TMPDIR=/scratch.global/coffm049/apptainer_build_tmp
export SINGULARITY_TMPDIR=/scratch.global/coffm049/apptainer_build_tmp
export APPTAINER_CACHEDIR=/scratch.global/coffm049/apptainer_build_tmp
export SINGULARITY_CACHEDIR=/scratch.global/coffm049/apptainer_build_tmp


apptainer build --fakeroot ancNreport.sif ancNreport.def
apptainer push ancNreport.sif oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest

apptainer build --fakeroot rfmix.sif rfmix.def
apptainer push rfmix.sif oras://ghcr.io/coffm049/gdcgenomicsqc/rfmix:latest

apptainer build --fakeroot mash.sif mash.def
apptainer push mash.sif oras://ghcr.io/coffm049/gdcgenomicsqc/mash:latest


apptainer build --fakeroot phenotypeSim.sif phenotypeSim.def
apptainer push phenotypeSim.sif oras://ghcr.io/coffm049/gdcgenomicsqc/phenotypesim:latest
```

# Module Load

## Structure (MSI format)

```
envs/gdcgenomicsqc/
└── 1.0         # TCL module definition (version file IS the module)
```

The module file must:
- Be named with the version number (e.g., `1.0`)
- Start with `#%Module`

## Testing

```bash
export MODULEPATH=/scratch.global/coffm049/GDCGenomicsQC/envs:$MODULEPATH
module avail gdcgenomicsqc
module load gdcgenomicsqc/1.0
module show gdcgenomicsqc/1.0
```


## What the module provides

When loaded via `module load gdcgenomicsqc/1.0`:

- Loads `apptainer` module (provides apptainer command)
- Adds `$basedir/bin` to PATH
- Sets `APPTAINER_CACHEDIR=/scratch.global/GDC/singularityimages`
- Sets `SNAKEMAKE_APPTAINER_PREFIX=/scratch.global/GDC/singularityimages`
- The `gdcgenomicsqc` wrapper runs snakemake with `--directory` set, so it works from any directory

These env vars allow snakemake to use cached apptainer images for offline execution.
