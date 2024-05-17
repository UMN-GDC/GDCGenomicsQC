<p align="center">
  <i>  A quality control pipeline for genomics data developed by the Masonic Institute of the Developing Brain at the University of Minnesota.</a></i>
  <br/>
</p>

![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54)
![R](https://img.shields.io/badge/r-%23276DC3.svg?style=for-the-badge&logo=r&logoColor=white)
![Shell Script](https://img.shields.io/badge/shell_script-%23121011.svg?style=for-the-badge&logo=gnu-bash&logoColor=white)

A quality control pipeline for genomics data developed by the Masonic Institute of the Developing Brain at the University of Minnesota. The pipeline is built utilizing Plink, Liftover, R-language, Python, and Bashed, and can be housed in a Docker image.

![GDC_pipeline_overview](https://github.com/UMN-GDC/GDCGenomicsQC/assets/140092486/e7f11909-9ab8-4def-90e5-c5f67c28a4bb)

## Standard Procedure *(Done in order)*

-   Exclude SNPs with greater than 10% missingness **(Plink)**

-   Exclude individuals with greater than 10% missingness **(Plink)**

-   Exclude SNPs with greater than 2% missingness **(Plink)**

-   Exclude individuals with greater than 2% missingness **(Plink)**

-   Compare sex assignments in input data set with imputed X chromosome coefficients **(Plink)**

    -   F-values \< 0.2 are assigned as female and F-values \> 0.8 are assigned as male others are flagged as problems and excluded from the dataset

-   Exclude SNPs with Minor Allele Frequency \< 0.1 **(Plink)**

-   Exclude SNPs where Hardy-Weinberg Equilibrium p-values \< 1e-6 for controls **(Plink)**

-   Exclude SNPs where Hardy-Weinberg Equilibrium p-values \< 1e-6 for cases **(Plink)**

-   Exclude SNPs that are highly coordinated using multiple correlation coefficients for a SNP regressed on all other SNPs simultaneously **(Plink)**

-   Exclude individuals with a parent-offspring relationship **(Plink)**

-   Exclude individuals with a pi_hat threshold \> 0.2 **(Primus)**

-   Principal Component Analysis (**FRAPOSA**)

## Optional QC Features

-   Synchronize data to GRCh38/hg38 (**Crossmap)**

-   Align strand orientation through (**Genotype Harmonizer)**

## Optional Pipeline Output

-   Gather information from log files created by standardized steps and create an automated report with tables and figures regarding each step.


# Installation
## Git clone
```shell
git clone https://github.com/coffm049/GDCGenomics/QC
```
### Reqirements
Still piecing these together.

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
