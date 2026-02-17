# building containers
Use the .SLRURM scripts

# Upload to GHCR
Register then push

```
apptainer registry login --username coffm049 docker://ghcr.io
apptainer push plink.sif oras://ghcr.io/coffm049/gdcgnomicsqc/plink:latest
```
