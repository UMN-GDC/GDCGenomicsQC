<p align="center">
  <i>  A quality control pipeline for genomics data developed by the Masonic Institute of the Developing Brain at the University of Minnesota.</a></i>
  <br/>
</p>

![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54)
![R](https://img.shields.io/badge/r-%23276DC3.svg?style=for-the-badge&logo=r&logoColor=white)
![Shell Script](https://img.shields.io/badge/shell_script-%23121011.svg?style=for-the-badge&logo=gnu-bash&logoColor=white)

A quality control pipeline for genomics data developed by the Masonic Institute of the Developing Brain at the University of Minnesota. The pipeline is built utilizing [Plink](https://www.cog-genomics.org/plink/), [Liftover](https://genome.ucsc.edu/cgi-bin/hgLiftOver), [R-language](https://www.r-project.org/), [Python](https://www.python.org/), and [bash](https://www.gnu.org/software/bash/), and  housed in a [Docker image](https://hub.docker.com/_/docker). The steps in the pipeline are detailed [here](https://gdcgenomicsqc.readthedocs.io/en/latest/)

## Quick Start

Choose the installation method that matches your environment:

<details>
<summary><b>HPC with Module (MSI/UMN) - Recommended</b></summary>

```bash
# Load the module
module use /projects/standard/gdc/public/envs/GDCGenomicsQC/envs
module load gdcgenomicsqc

# Set up Snakemake (if not already available)
conda config --add envs_dirs /projects/standard/gdc/public/envs
conda activate snakemake

# Clone and run
git clone https://github.com/UMN-GDC/GDCGenomicsQC.git
cd GDCGenomicsQC/workflow
snakemake --profile ../profiles/sandbox --configfile /path/to/your/config.yaml
```
</details>

<details>
<summary><b>HPC without Module</b></summary>

```bash
# Clone the repository
git clone https://github.com/UMN-GDC/GDCGenomicsQC.git
cd GDCGenomicsQC

# Set up Snakemake environment
conda env create -n snakemake -f envs/snakemake.yml
conda activate snakemake

# Run with HPC profile
cd workflow
snakemake --profile ../profiles/hpc --configfile /path/to/your/config.yaml
```
</details>

<details>
<summary><b>Interactive/Local</b></summary>

```bash
# Clone the repository
git clone https://github.com/UMN-GDC/GDCGenomicsQC.git
cd GDCGenomicsQC

# Set up Snakemake environment
conda env create -n snakemake -f envs/snakemake.yml
conda activate snakemake

# Run interactively
cd workflow
snakemake --cores=4 --use-conda --configfile /path/to/your/config.yaml
```
</details>

<details>
<summary><b>Singularity/Apptainer Only</b></summary>

```bash
# Install Snakemake via pip
pip install snakemake snakemake-executor-plugin-slurm

# Clone the repository
git clone https://github.com/UMN-GDC/GDCGenomicsQC.git
cd GDCGenomicsQC

# Run with Singularity
cd workflow
snakemake --use-singularity --profile ../profiles/hpc --configfile /path/to/your/config.yaml
```
</details>

---

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
- Validate your configuration against the schema:
  ```bash
  snakemake --validate --configfile config/config.yaml
  ```
- Snakemake expects you to execute from `GDCGenomicsQC/workflow` 
- Run the desired workflow (by default looks in `config/config.yaml`) using the `--configfile=<path/to/confi.yaml/` flag

To have SLURM dispatch it without dependency on your terminal being open these snakemake calls can be called in a SLURM script.
 An example is stored at workflow/example.SLURM

### Running on a New HPC (No Module Available)

If you're on an HPC system that doesn't have the GDCGenomicsQC module installed, you'll need to set it up manually:

#### 1. Clone the Repository
```bash
git clone https://github.com/UMN-GDC/GDCGenomicsQC.git
cd GDCGenomicsQC
```

#### 2. Set Up Snakemake Environment
```bash
# Option A: Create a conda/mamba environment
conda env create -n snakemake snakemake snakemake-executor-plugin-slurm
conda activate snakemake

# Option B: Use existing snakemake if available
module load snakemake
```

#### 3. Configure and Run

**Option A: Standard HPC profile (requires internet):**
```bash
cd GDCGenomicsQC/workflow
snakemake --profile ../profiles/hpc --configfile /path/to/your/config.yaml
```

**Option B: Sandbox profile (if pre-cached images available):**
```bash
cd GDCGenomicsQC/workflow
snakemake --profile ../profiles/sandbox --configfile /path/to/your/config.yaml
```

**Option C: Interactive profile (for testing):**
```bash
cd GDCGenomicsQC/workflow
snakemake --profile ../profiles/interactive --configfile /path/to/your/config.yaml
```

#### Requesting Module Installation

If you'd like the module installed system-wide, contact your HPC administrators with:
- The path to the repository: `/path/to/GDCGenomicsQC`
- The module location: `/path/to/GDCGenomicsQC/envs/gdcgenomicsqc`
- The wrapper script: `/path/to/GDCGenomicsQC/envs/bin/gdcgenomicsqc`

## Detailed usage 
After cloning this repository the steps to run this pipeline are as follows:
1.	To run pipeline with SLURM scheduler (reccomended): `snakemake --profile=../profiles/hpc`
2.	To run pipeline interactively: `snakemake --profile=../profiles/interactive`

These profiles specify using the singularity images (`--use-singularity`), but if desired you can run them with the `--use-conda` flag, which will construct and cache the conda envs locally. Just note that this does not work well with SLURM schedulers (`--executor slurm`), but will work fine when running interactively.

### Sandbox Profile (Offline/Isolated Environment)

The `sandbox` profile is designed for environments with limited or no internet access. It uses a pre-cached set of Singularity images stored in a shared location.

#### Installation

1. **First-time setup** (requires internet to cache images):
   ```bash
   cd GDCGenomicsQC/workflow
   snakemake --profile ../profiles/sandbox -n
   ```

   This will pull all required Singularity images to the cache at `/scratch.global/GDC/singularityimages/`.

2. **Configure module** (optional, for easier access):
   ```bash
   # If using the module file from the repo:
   module use /path/to/GDCGenomicsQC/envs
   module load gdcgenomicsqc
   ```

#### Running with Sandbox Profile

**Using the module (recommended):**
```bash
module use /path/to/GDCGenomicsQC/envs
module load gdcgenomicsqc

# Verify environment is set up
echo $SINGULARITY_CACHEDIR
echo $SNAKEMAKE_SINGULARITY_PREFIX

# Run from any directory - --directory is set automatically
gdcgenomicsqc --configfile /path/to/your/config.yaml
```

**Without module:**
```bash
cd GDCGenomicsQC/workflow
snakemake --profile ../profiles/sandbox --configfile /path/to/your/config.yaml
```

#### Comparison: Sandbox vs HPC/Interactive vs New HPC

| Feature | Sandbox | HPC | New HPC (No Module) | Interactive |
|---------|---------|-----|---------------------|-------------|
| Container cache | Pre-cached at `/scratch.global/GDC/singularityimages` | Downloads on demand | Downloads on demand | Downloads on demand |
| Internet required | No (after initial setup) | Yes | Yes | Yes |
| Profile location | `profiles/sandbox` | `profiles/hpc` | `profiles/hpc` or `profiles/sandbox` | `profiles/interactive` |
| Use case | Offline/air-gapped systems | Standard HPC runs | Fresh HPC without module | Local testing |
| Module available | Yes (`gdcgenomicsqc`) | Yes (system) | No - manual setup required | No |
| Setup required | Clone repo, run once to cache | Clone repo | Clone + create conda env | Clone + create conda env |

**Setup time comparison:**
- **Sandbox**: ~5 min initial setup + first run to cache (internet required), then offline
- **HPC (with module)**: ~2 min to load module
- **New HPC (no module)**: ~10-15 min (clone + conda setup)
- **Interactive**: ~10-15 min (clone + conda setup)

**Switching between profiles:**
```bash
# Use sandbox (offline)
snakemake --profile ../profiles/sandbox ...

# Use HPC (standard)
snakemake --profile ../profiles/hpc ...

# Use interactive (local)
snakemake --profile ../profiles/interactive ...
```
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
Details for each specific project are configured in a .yaml file. An example is provided in `GDCGenomicsQC/config/config.yaml`

### Core Parameters
```yaml
INPUT_FILE: "/path/to/input/vcf"
OUT_DIR: "/path/to/output"
REF: "/path/to/reference"
vcf_template: "/path/to/vcf_template_chr{chr}.vcf.gz"
local-storage-prefix: "/path/to/.snakemake/storage"
```

### QC Options
```yaml
relatedness:
    method: "0"  # 0=assume unrelated, 1=KING, 2=PRIMUS

SEX_CHECK: false
GRM: true  # Generate genetic relationship matrix

thin: true  # Thin input data
```

### Ancestry Estimation
```yaml
ancestry:
    threshold: 0.8  # Minimum ancestry proportion
    model: "pca"  # Options: pca, umap, vae, rfmix

internalPCA:
    plot: true
    color_by: null  # Column to color by (e.g., "ancestry", "sex")
    phenotype_file: null
```

### Local Ancestry (RFMIX)
```yaml
localAncestry:
  RFMIX: true
  test: true
  thin_subjects: 0.1
  figures: "figures"
```

### Phenotype Simulation
```yaml
phenotypeSimulation:
    ancestries: ["AFR", "EUR"]
    n_sims: 10
    heritability: 0.4
    rho: 0.8
    maf: 0.05
    seed: 42
    skip_thinning: true
    thin_count_snps: 1000000
    thin_count_inds: 10000
```

### SNP Heritability
```yaml
snpHerit:
    pheno: null  # Path to phenotype file (required)
    covar: null  # Path to covariate file(s), comma-separated for multiple files (required)
    out: phenoEsts
    method: "AdjHE"  # AdjHE, GCTA, PredLMM, SWD, COMBAT
    random_groups: false
    npc: [10]
    loop_covs: false
    mpheno: 1
```

### VCF Conversion
```yaml
convertNfilt:
    info_r2_min: null  # Minimum R2 from INFO field
    filter_pass: true  # Only keep FILTER==PASS
    qual_min: null     # Minimum QUAL value
```

### Environment
```yaml
conda-frontend: mamba
```

## Available Rules
The pipeline includes the following rules:

| Rule | Description |
|------|-------------|
| `convertNfilt` | Convert VCF to PLINK format with filtering |
| `Initial_QC` | Initial quality control (MAF, missingness) |
| `Standard_QC` | Standard QC (HWE, inversion regions, sex check) |
| `Relatedness` | Check and filter related samples |
| `PCAreference` | PCA on reference panel and sample projection |
| `UMAP` | UMAP dimensionality reduction |
| `estimateAncestry` | Estimate global ancestry using PC/UMAP |
| `classifyAncestry` | Classify ancestry with multiple models |
| `Phase` | Phase genotypes with SHAPEIT4 |
| `RFMIX` | Local ancestry inference with RFMIX |
| `simulatePhenotype` | Simulate phenotypes for GWAS |
| `snpHerit` | Estimate SNP heritability |

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
    - ref.<eigenvec, eigenval> - PCA information on the reference panel
    - <ref, sample>RefPCscores.sscore - projection of sample and reference on the reference PCs
    - umap_<sample, ref>.csv - the UMAP embeddings of the sample and reference
    - posterior_probabilities.tsv - posterior probability for each ancestry
    - sample_coords.tsv - sample coordinates in latent space
    - ref_coords.tsv - reference coordinates in latent space
    - posterior_probability_stacked_<model>.svg - ridge plot of ancestry proportions
 - 02-localAncestry
    - chr<chr>.lai.<fb, msp, rfmix.Q, sis> - rfmix output specifying posterior probability, most probable LAI label, and chromosome aggregated admixing proportion
    - chr<chr>phased.vcf.gz - the Shapeit4 phased haplotypes
 - 03-snpHeritability - SNP heritability estimates

## Create DAG
```bash
# snakemake run_pca --dag --configfile ../config/abcd_config.yaml 2>/dev/null | dot -Tpng > /scratch.global/coffm049/GDCGenomicsQC/dag_run_pca.png
# snakemake run_pca --rulegraph --configfile ../config/abcd_config.yaml 2>/dev/null | dot -Tpng > /scratch.global/coffm049/GDCGenomicsQC/dag_run_pca.png
snakemake run_pca --dag mermaid-js --configfile ../config/abcd_config.yaml 2>/dev/null > ../docs/dag_pca.mmd
snakemake run_pca --rulegraph mermaid-js --configfile ../config/abcd_config.yaml 2>/dev/null > ../docs/rulegraph_pca.mmd
snakemake run_initialQC --dag mermaid-js --configfile ../config/abcd_config.yaml 2>/dev/null > ../docs/dag_initialQC.mmd
snakemake run_initialQC --rulegraph mermaid-js --configfile ../config/abcd_config.yaml 2>/dev/null > ../docs/rulegraph_initialQC.mmd

```

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
