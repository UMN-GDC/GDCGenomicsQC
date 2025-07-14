<p align="center">
  <i>  A quality control pipeline for genomics data developed by the Masonic Institute of the Developing Brain at the University of Minnesota.</a></i>
  <br/>
</p>

![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54)
![R](https://img.shields.io/badge/r-%23276DC3.svg?style=for-the-badge&logo=r&logoColor=white)
![Shell Script](https://img.shields.io/badge/shell_script-%23121011.svg?style=for-the-badge&logo=gnu-bash&logoColor=white)

A quality control pipeline for genomics data developed by the Masonic Institute of the Developing Brain at the University of Minnesota. The pipeline is built utilizing Plink, Liftover, R-language, Python, and Bashed, and can be housed in a Docker image.

# Installation
## Git clone
```shell
git clone https://github.com/UMN-GDC/GDCGenomicsQC.git
```
### Reqirements
-	Access to MSI computing resources.
-	Genomic files in plink bed formatting (bim, bed, & fam)

# Usage
After cloning this repository the steps to run this pipeline are as follows:
1.	To run pipeline: `sh ./GDCGenomicsQC/Run.sh`
2.	Flags to be appended to run command `sh ./GDCGenomicsQC/Run.sh --flag1 option1 --flag2 option2`
 -	`--set_working_directory`	Provide a path to where you'd like the outputs to be stored
 -	`--input_directory`	Provide the path to where the bim/bed/fam data is stored
 -	`--input_file_name`	Provide a path to where you'd like the outputs to be stored
 -	`--path_to_github_repo`	Provide the path to where the bim/bed/fam data is stored
 -	`--user_x500`	Provide a path to where you'd like the outputs to be stored
 -	`--use_crossmap`	Provide the path to where the bim/bed/fam data is stored
 -	`--use_genome_harmonizer`	Provide a path to where you'd like the outputs to be stored
 -	`--use_genome_harmonizer`	Provide the path to where the bim/bed/fam data is stored
 -	`--use_genome_harmonizer`	Provide the path to where the bim/bed/fam data is stored
 -	`--use_genome_harmonizer`	Provide the path to where the bim/bed/fam data is stored
 -	`--use_genome_harmonizer`	Provide the path to where the bim/bed/fam data is stored
3.	Execute or source the copy of the settings_file_template.sh to create a file ending with _wrapper.sh.
4.	Run the wrapper file created as an sbatch submission. Below is an example of how to do so. 
	```shell
	sbatch sample_wrapper.sh
	```

![GDC_pipeline_overview](https://github.com/UMN-GDC/GDCGenomicsQC/assets/140092486/e7f11909-9ab8-4def-90e5-c5f67c28a4bb)

## Standard Procedure *(Done in order)*

### Module 1: Crossmap (optional)

This moudle converts the genome build to GRCh38.  The default setting is conversion from GRCh37 to GRCh38.  This is controlled using the `--use_crossmap` flag.  The default setting is `1`, but if the build is already GRCh38 this flag should be manually set to `0`.

### Module 2: GenotypeHarmonizer (optional)

This module aligns sample to the GRCh38 reference genome.  The reference genome alignment file we use is `ALL.hgdp1kg.filtered.SNV_INDEL.38.phased.shapeit5.vcf` which is designed for the GRCh38 build.  All subsequent steps in this pipeline assume we have an aligned GRCH38 build study sample.  We can control for alignment using the flag `--use_genome_harmonizer`.  The defualt flag setting is `1`, but this can be manually set to `0`.

### Module 3: Initial QC

We perform an initial round of Quality Control prior to runing relatedness checks.

-   Exclude SNPs with greater than 10% missingness **(Plink)**
-   Exclude individuals with greater than 10% missingness **(Plink)**
-   Exclude SNPs with greater than 2% missingness **(Plink)**
-   Exclude individuals with greater than 2% missingness **(Plink)**

### Module 4: Relatedness

We first run kinship test using KING.  This separates the plink dataset into related and unrelated plink study samples.  Many of the subsequent steps involving ancestry estimation require the sample subjects to be unrelated.  To provide some analysis of related subjects we also perform PC-AiR and PC-Relate in this module.  We convert PLINK format to GDS format and perform an initial kinship estimate with the KING algorithm, but pivot to perform Principal Component Analysis and infer genetic ancestry on each inidvidual.  We finally run PC-Relate to calculate highly accurate measures of genetic relatedness. (kinship coefficients, IBD probabilities).

### Module 5: Standard QC

This runs the standard quality control measures expected of GWAS on unrelated individuals.  To compensate for potential related individuals, we perform QC measures on the unrelated study samples from Module 4: Relatedness, but apply similar filtering standards on the related individuals.  We create the unrelated QC dataset, but also a related QC dataset which has the same SNPs extracted as the unrleated dataset.  The following are the standard QC steps.  Custom QC steps can be provided by generating a custom script.

-   Exclude SNPs with greater than 10% missingness **(Plink)**
-   Exclude individuals with greater than 10% missingness **(Plink)**
-   Exclude SNPs with greater than 2% missingness **(Plink)**
-   Exclude individuals with greater than 2% missingness **(Plink)**
-   Compare sex assignments in input data set with imputed X chromosome coefficients **(Plink)**
    -   F-values \< 0.2 are assigned as female and F-values \> 0.8 are assigned as male others are flagged as problems and excluded from the dataset
-   Exclude SNPs with Minor Allele Frequency \< 0.01 **(Plink)**
-   Exclude SNPs where Hardy-Weinberg Equilibrium p-values \< 1e-6 for controls **(Plink)**
-   Exclude SNPs where Hardy-Weinberg Equilibrium p-values \< 1e-10 for cases **(Plink)**
-   Exclude SNPs that are highly coordinated using multiple correlation coefficients for a SNP regressed on all other SNPs simultaneously **(Plink)**
-   Exclude individuals with a parent-offspring relationship **(Plink)**
-   Exclude individuals with a pi_hat threshold \> 0.2 **(Primus)**
-   Principal Component Analysis (**FRAPOSA**)

### Module 6: Phasing

We perform phasing using shapeit4.2.  This is necessary to use rfmix to infer ancestry.  This module first recodes the QC unrleated dataset into vcf format separated by chromosome.  Phasing is then performed using reference map `chr${CHR}.b38.gmap.gz`.  The files are then ready to be run using rfmix.

### Module 7: rfmix

We infer ancestry of individual samples using rfmix.  In addition to phased file provide by `module 6: phasing`, we also need reference genome that has also been phased, a population map or super population map file, and a genetic map file for GRch 38 build.  In our script we use reference genome `hg38_phased.vcf.gz`, super population map file `super_population_map_file.txt`, and genetic map `genetic_map_hg38.txt`.  Initially we generate posterior probabilities for each ancestry by sample.  These posterior probabilities represent each of the 22 chromosomes.  To get global ancestries for each individual, we take the mean posterior probabilities for each super population across all 22 .Q files.  Assignment to globabl ancestry is based on the highest posterior probabilty that is greater than 0.8.  If no posterior probablitiy is greater than 0.8, that subject classified as `Other`.

### Module 8: population stratification

We separate the samples into individual plink files based on their most probable posterior ancestries dtermined in Module 7: rfmix.

### Module 9: ancestry plots

This module provides visualization for ancestry estimation.  We provide two sets of plots.  

-   GAP:  This visualization the proportion of each ancestry in individual samples.
-   LAP:  This visualization shows most probable posterior ancestry by regions of the chromosome.

### Subpopulation QC

We can also provide individual QC steps stratified by population.

## Optional Pipeline Output

-   Gather information from log files created by standardized steps and create an automated report with tables and figures regarding each step.



## Docker ![Docker](https://img.shields.io/badge/docker-%230db7ed.svg?style=for-the-badge&logo=docker&logoColor=white)
Still in development


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
