<p align="center">
  <i>  A quality control pipeline for genomics data developed by the Masonic Institute of the Developing Brain at the University of Minnesota.</a></i>
  <br/>
</p>

![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54)
![R](https://img.shields.io/badge/r-%23276DC3.svg?style=for-the-badge&logo=r&logoColor=white)
![Shell Script](https://img.shields.io/badge/shell_script-%23121011.svg?style=for-the-badge&logo=gnu-bash&logoColor=white)

A quality control pipeline for genomics data developed by the Masonic Institute of the Developing Brain at the University of Minnesota. The pipeline is built utilizing [Plink](https://www.cog-genomics.org/plink/), [Liftover](https://genome.ucsc.edu/cgi-bin/hgLiftOver), [R-language](https://www.r-project.org/), [Python](https://www.python.org/), and [bash](https://www.gnu.org/software/bash/), and  housed in a [Docker image](https://hub.docker.com/_/docker). The steps in the pipeline are detailed [here](https://gdcgenomicsqc.readthedocs.io/en/latest/)



## Features
- State-of-the-art genomics quality control pipeline
    - Assesses relatedness
    - Assess global and local ancestry
    - Controls for relatedness and genetic ancestry in QC steps
    - SNP-heritability methods for multiple ancestries
    - PRS methods for multiple ancestries
    - Easy extensibility, reproducibility, and modularity

- Workflow management with [Snakemake](https://snakemake.github.io/)
    - Smart execution of workflow steps
        - Specify desired output and Snakemake will back-construct necessary steps to create it
    - Controlled conda environments automatically handled 
    - Controlled containers automatically handled (under construction)
    - Workflow handling on local computers and with SLURM scheduling
    - Automated report generation


# Usage 
Requirements:
- Access to HPC computing resources with SLURM scheduler (though it can still run in any terminal , just --executor slurm won't function).
- Snakemake
    - Can be installed with `conda env create -n snakemake snakemake snakemake-executor-plugin-slurm conda`
	- If you are running on MSI at UMN this environment already exists and you won't need to reinstall it 
	- You can simply add the list of GDC conda envs by running `conda config --add envs_dirs /projects/standard/gdc/public/envs`
    - this installs the conda environment called snakemake
    - Activate conda env: `conda activate snakemake`

## Installation

```shell
git clone https://github.com/UMN-GDC/GDCGenomicsQC.git
cd GDCGenomicsQC
conda env create -n snakemake snakemake snakemake-executor-plugin-slurm
```

## Using Snakemake workflows
- update config files as necessary (located at `config/config.yaml`)
    - Update Inputs and outputs and methods details as desired
    - updated SLURM group name for accounting purposes
- Snakemake expects you to execute from `GDCGenomicsQC/workflow` 
- Run the desired workflow (by default looks in `config/config.yaml`) using the `--configfile=<path/to/confi.yaml/` flag

To have SLURM dispatch it without dependency on your terminal being open these snakemake calls can be called in a SLURM script.
 An example is stored at workflow/example.SLURM

## Detailed usage 
After cloning this repository the steps to run this pipeline are as follows:
1.	To run pipeline with SLURM scheduler (reccomended): `snakemake --profile=../profiles/hpc`
2.	To run pipeline interactively: `snakemake --profile=../profiles/interactive`

These profiles specify using the singularity images (`--use-singularity`), but if desired you can run them with the `--use-conda` flag, which will construct and cache the conda envs locally. Just note that this does not work well with SLURM schedulers (`--executor slurm`), but will work fine when running interactively.
 - `--configfile </path/to/configfile>` path to .yaml configuring your desired run
 - to execute it somewhere else add these flags `--directory /path/to/GDCGenomicsQC/workflow --snakefile /path/to/GDCGenomicsQC/workflow/Snakefile`
 - For older versions of snakemake (if you dindn't install conda env create snakame as specified above) run with a slurm scheduler by appending `--cluster "sbatch --parsable"`
 - `--jobs` maximum number of slurm jobs to submit at once. If this is smaller that 22, note that steps that run per autosomal chromosome will be submitted sequentially in phases
 - `<Rule Name>` if you only want to run specific aspects of the pipeline you can specify the rule you want to run through
    - Initial_QC
    - PCA
    - RFMIX
 - `--report --report-stylesheet /path/to/GDCGenomicsQC/report/stylesheet.css` tells snakemake to create a summary .html report at `GDCGenomicsQC/workflow/report.html`

### Recommended calls

An example batch job is included at `workflow/example.SLURM` for easier adaptation to your workflow. If available, we reccomend letting SLURM handle the disbatching and generating a report.

```bash
snakemake --profile=../profiles/hpc \
    --report --report-stylesheet /path/to/GDCGenomicsQC/report/stylesheet.css \
    --configfile </path/to/config.yaml> --directory </path/to/GDCGenomicsQC/workflow> --snakefile </path/to/GDCGenomicsQC/workflow/Snakefile>
```

For local execution
```bash
snakemake --cores=4 --use-conda \
    --configfile </path/to/config.yaml> --directory </path/to/GDCGenomicsQC/workflow> --snakefile </path/to/GDCGenomicsQC/workflow/Snakefile>
```


For running upt a to a certain point (i.e. PCA)
```bash
snakemake --cores=4 --use-conda \
    --configfile </path/to/config.yaml> --directory </path/to/GDCGenomicsQC/workflow> --snakefile </path/to/GDCGenomicsQC/workflow/Snakefile> \
    PCA
```

As of recent Snakemake bug report this has been incorporated into the `hpc` profile and might be removed in future versions
```bash
snakemake --executor=slurm --use-singularity --local-storage $(pwd)/.snakemake/storage
```


## Configuration
Details for each specific projcet are configured in a .yaml file. An example is provided in `GDCGenomicsQC/config/config.yaml`
```yaml
INPUT_FILE: "/projects/standard/gdc/public/Ref/toyData/1kgSynthetic"
OUT_DIR: "/scratch.global/coffm049/toyPipeline"
REF: "/projects/standard/gdc/public/Ref"

# Tool-specific parameters
relatedness:
    method: "0"

SEX_CHECK : false
GRM : false
RFMIX : true
rfmix_test : true

thin: true
```

![GDC_pipeline_overview](https://github.com/UMN-GDC/GDCGenomicsQC/assets/140092486/e7f11909-9ab8-4def-90e5-c5f67c28a4bb)

# Output
The output directory is organized as follows
- Genomic data derivatives are prefixed with the number in which they are ran
	- So far this includes 01-globalAncestry, 02-localAncestry
- Unnumbered directories include
	- simulations - simulated phenotypes where the subdirectories describe the combination of ancestires included in each given simulation
	- data subsets such as `full` and identified ancestry subsets. for thousand genome reference panel this includes `AFR`, `AMR`, `EAS`, `EUR`, `SAS`
    - Each of these direcotries contain QC output assumning each is a homogeneous group.
    	- initialFilter_<chr> - are filtered for MAF, and variant missingness
        - initialFilter - Are fully combined genomes additionally with sample missingness filter
        - standardFilter - additionally filters for inversion regions, hardy-weinberg equilibrium, and check's sex (if specified in the config)
        - standardFilter.LDpruned - additionally is filtered for linkage diseqilibrium and the specified level (default is 500 10 0.1)
 - 01-globalAncestry
 	- ref.<acount, eigenval, eigenvec, eigenvec.allele> - PCA information on the reference panel
  	- <ref, sample>RefPCscores.sscore - projection of sample and reference  onthe reference PCs
    - umap_<sample, ref>.csv - the UMAP embeddings of the sample and refrence
    - latentDistantRelatedness.tsv - combination of PC's, UMAPS, and most probable label based on random forest training on PC and UMAP embeddings respectively
    - RF<pc, umap>.Rds - R object containing the trained random forest on specified latent variables
- 02-localAncestry
  	- chr<chr>.lai.<fb, msp, rfmix.Q, sis> - rfmix output specifying posterior probaiblity, most probable LAI label, and chromosome aggregated admixing proportion
  	- chr<chr>phased.vcf.gz - the Shapeit4 phased haplotypes
## Contributing

GDCGenomicsQC is built and maintained by a small team – we'd love your help to fix bugs and add features!

Before submitting a pull request _please_ discuss with the core team by creating or commenting in an issue on [GitHub](https://www.github.com/coffm049/GDCGenomics/issues) – we'd also love to hear from you in the [discussions](https://www.github.com/coffm049/GDCGenomics/discussions). This way we can ensure that an approach is agreed on before code is written. This will result in a much higher likelihood of your code being accepted.

If you’re looking for ways to get started, here's a list of ways to help us improve:

- Issues with [`good first issue`](https://github.com/outline/outline/labels/good%20first%20issue) label
- Developer happiness and documentation
- Bugs and other issues listed on GitHub

## Tests
This is still under construction


# License

[MIT licensed](LICENSE).
