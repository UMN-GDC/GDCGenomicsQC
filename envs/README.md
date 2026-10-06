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


apptainer build --fakeroot genomeUtils.sif genomeUtils.def
apptainer push genomeUtils.sif oras://ghcr.io/coffm049/gdcgenomicsqc/genomeutils:latest

apptainer build --fakeroot rfmix.sif rfmix.def
apptainer push rfmix.sif oras://ghcr.io/coffm049/gdcgenomicsqc/rfmix:latest

apptainer build --fakeroot mash.sif mash.def
apptainer push mash.sif oras://ghcr.io/coffm049/gdcgenomicsqc/mash:latest


apptainer build --fakeroot phenotypeSim.sif phenotypeSim.def
apptainer push phenotypeSim.sif oras://ghcr.io/coffm049/gdcgenomicsqc/phenotypesim:latest

apptainer build --fakeroot popvae.sif popvae.def
apptainer push popvae.sif oras://ghcr.io/coffm049/gdcgenomicsqc/popvae:latest

apptainer build --fakeroot --build-arg TOKEN_FILE=~/.config/gdcgenomicsqc/predlmmace.token predlmmAce.sif predlmmAce.def
apptainer push predlmmAce.sif oras://ghcr.io/coffm049/gdcgenomicsqc/predlmmace:latest
```

## predlmmAce: private dependency

`envs/predlmmAce.yml` installs predlmm-ace from `git+https://github.com/SAONLIB/predlmm-ace.git`, and that repo is private (anonymous clone returns 404). The Apptainer build therefore needs a GitHub token with read access to it.

Create the token file once, outside the repo:

```bash
mkdir -p ~/.config/gdcgenomicsqc
printf '%s' "$GITHUB_TOKEN" > ~/.config/gdcgenomicsqc/predlmmace.token
chmod 600 ~/.config/gdcgenomicsqc/predlmmace.token
```

Then build with `--build-arg TOKEN_FILE=<path>`, as shown above.

**Pass the path, never the token.** Anything in `%arguments` is substituted into `%post` and stored verbatim in the SIF's embedded definition, where `apptainer inspect --deffile predlmmAce.sif` will print it. `--build-arg TOKEN=<secret>` would publish the token with the image. Build-time environment variables are not an alternative either: they do not reach `%post`.

`predlmmAce.def` reads the token from `TOKEN_FILE` and hands it to git through an ephemeral `GIT_CONFIG_*` credential helper, so the secret is never written into the image, never appears in the def or in `ps`, and is unset before the filesystem is squashed. Verified against the built image: the token does not appear in the embedded definition, in the SIF bytes, or anywhere in the extracted filesystem, and the image contains no `.gitconfig`, `.git-credentials`, or `.netrc`.

### Published image

`oras://ghcr.io/coffm049/gdcgenomicsqc/predlmmace:latest`

- digest: `sha256:3b31bedfe4926f0c4cbb3803fe99e5c125c49e4dd8da6ab2683961f4121f5d60`
- size: 165MB
- provides `predlmm-fit`, `predlmm-profile-se`, `predlmm-grm-to-nystrom`, `predlmm-select-knots`, and `gcta` v1.94.1

Round-tripped from GHCR: `apptainer pull` returns the same digest, and `workflow/scripts/run_predlmm_ace.py --help` runs inside the pulled image. `workflow/rules/snpHeritRelated.smk` references this URI in a `container:` directive alongside the existing `conda:` directive, so the same rules work with or without `--software-deployment-method apptainer`.

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
