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
```
