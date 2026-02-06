<p align="center">
  <i>  A quality control pipeline for genomics data developed by the Masonic Institute of the Developing Brain at the University of Minnesota.</a></i>
  <br/>
</p>

![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54)
![R](https://img.shields.io/badge/r-%23276DC3.svg?style=for-the-badge&logo=r&logoColor=white)
![Shell Script](https://img.shields.io/badge/shell_script-%23121011.svg?style=for-the-badge&logo=gnu-bash&logoColor=white)

A quality control pipeline for genomics data developed by the Masonic Institute of the Developing Brain at the University of Minnesota. The pipeline is built utilizing [Plink](https://www.cog-genomics.org/plink/), [Liftover](https://genome.ucsc.edu/cgi-bin/hgLiftOver), [R-language](https://www.r-project.org/), [Python](https://www.python.org/), and [bash](https://www.gnu.org/software/bash/), and  housed in a [Docker image](https://hub.docker.com/_/docker).

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
conda env create -n snakemake snakemake snakemake-executor-plugin-slurm conda
```

## Using Snakemake workflows
- update config files as necessary (located at config/config.yaml)
- Run the desired workflow (by default looks in config/config.yaml)
- Snakemake expects you to execute from `GDCGenomicsQC/workflow` 

To have SLURM dispatch it without dependency on your terminal being open these snakemake calls can be called in a SLURM script.
 An example is stored at workflow/example.SLURM

## Detailed usage 
After cloning this repository the steps to run this pipeline are as follows:
1.	To run pipeline: `snakemake --cores=4`
    - `--cores` (required) specify maximum number of threads for each step. Each step will use up to this number of threads, but each has their own internally specified number of threads to execute
2.	Flags to be appended to run command `snakemake --cores=4`
 - `--use-conda` has snakemake construct the conda envs  in `GDCGenomicsQC/envs` and cache them for running
 - `--use-singularity` has snakemake construct the conda envs  in `GDCGenomicsQC/envs` and cache them for running
 - `--configfile </path/to/configfile>` path to .yaml configuring your desired run
 - `--executor slurm`
 - to execute it somewhere else add these flags `--directory /path/to/GDCGenomicsQC/workflow --snakefile /path/to/GDCGenomicsQC/workflow/Snakefile`
    - For older versions of snakemake (if you dindn't install conda env create snakame as specified above) this is `--cluster "sbatch --parsable"`
 - `--jobs` maximum number of slurm jobs to submit at once. If this is smaller that 22, note that steps that run per autosomal chromosome will be submitted sequentially in phases
 - `<Rule Name>` if you only want to run specific aspects of the pipeline you can specify the rule you want to run through
    - Initial_QC
    - PCA
    - RFMIX
 - `--report --report-stylesheet /path/to/GDCGenomicsQC/report/stylesheet.css` tells snakemake to create a summary .html report at `GDCGenomicsQC/workflow/report.html`

### Recommended calls

An example batch job is included at `workflow/example.SLURM` for easier adaptation to your workflow. If available, we reccomend letting SLURM handle the disbatching and generating a report.
```bash
snakemake --cores=4 --use-conda \
    --configfile </path/to/config.yaml> --directory </path/to/GDCGenomicsQC/workflow> --snakefile </path/to/GDCGenomicsQC/workflow/Snakefile> \
    --executor slurm \
    --jobs 25 \
    --report --report-stylesheet /path/to/GDCGenomicsQC/report/stylesheet.css
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

For generating the report
```bash
snakemake --cores=4 --use-conda \
    --configfile </path/to/config.yaml> --directory </path/to/GDCGenomicsQC/workflow> --snakefile </path/to/GDCGenomicsQC/workflow/Snakefile> \
    --report --report-stylesheet /path/to/GDCGenomicsQC/report/stylesheet.css


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



# Standard Procedure *(Done in order)*

Full documentation is [here](https://docs.google.com/document/d/147DG4Po_TaSFyiXkmRlszrmQIjsp4UNNe_u90slZbD8/edit?tab=t.0). Each step is contained in `workflow/rules`

### Module 1: Crossmap (optional)

This module converts the genome build to GRCh38.  The default setting is conversion from GRCh37 to GRCh38.  This is controlled using the `--use_crossmap` flag.  The default setting is `1`, but if the build is already GRCh38 this flag should be manually set to `0`.

### Module 2: GenotypeHarmonizer (optional)

This module aligns sample to the GRCh38 reference genome.  The reference genome alignment file we use is `ALL.hgdp1kg.filtered.SNV_INDEL.38.phased.shapeit5.vcf` which is designed for the GRCh38 build.  All subsequent steps in this pipeline assume we have an aligned GRCH38 build study sample.  We can control for alignment using the flag `--use_genome_harmonizer`.  The default flag setting is `1`, but this can be manually set to `0`.

### Module 3: Initial QC

We perform an initial round of Quality Control prior to running relatedness checks.

-   Exclude SNPs with greater than 10% missingness **(Plink)**
-   Exclude individuals with greater than 10% missingness **(Plink)**
-   Exclude SNPs with greater than 2% missingness **(Plink)**
-   Exclude individuals with greater than 2% missingness **(Plink)**

### Module 4: Relatedness

We first run kinship test using KING.  This separates the Plink dataset into related and unrelated Plink study samples.  Many of the subsequent steps involving ancestry estimation require the sample subjects to be unrelated.  To provide some analysis of related subjects we also perform PC-AiR and PC-Relate in this module.  We convert PLINK format to GDS format and perform an initial kinship estimate with the KING algorithm, but pivot to perform Principal Component Analysis and infer genetic ancestry on each individual.  We finally run PC-Relate to calculate highly accurate measures of genetic relatedness. (kinship coefficients, IBD probabilities).

### Module 5: Standard QC

This runs the standard quality control measures expected of GWAS on unrelated individuals.  To compensate for potential related individuals, we perform QC measures on the unrelated study samples from Module 4: Relatedness, but apply similar filtering standards on the related individuals.  We create the unrelated QC dataset, but also a related QC dataset which has the same SNPs extracted as the unrelated dataset.  The following are the standard QC steps.  Custom QC steps can be provided by generating a custom script.

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

We perform phasing using shapeit4.2.  This is necessary to use rfmix to infer ancestry.  This module first recodes the QC unrelated dataset into vcf format separated by chromosome.  Phasing is then performed using reference map `chr${CHR}.b38.gmap.gz`.  The files are then ready to be run using rfmix.

### Module 7: rfmix

We infer ancestry of individual samples using rfmix.  In addition to phased file provide by `module 6: phasing`, we also need reference genome that has also been phased, a population map or super population map file, and a genetic map file for GRch 38 build.  In our script we use reference genome `hg38_phased.vcf.gz`, super population map file `super_population_map_file.txt`, and genetic map `genetic_map_hg38.txt`.  Initially we generate posterior probabilities for each ancestry by sample.  These posterior probabilities represent each of the 22 chromosomes.  To get global ancestries for each individual, we take the mean posterior probabilities for each super population across all 22 .Q files.  Assignment to global ancestry is based on the highest posterior probability that is greater than 0.8.  If no posterior probability is greater than 0.8, that subject classified as `Other`.

### Module 8: ancestry plots

This module provides visualization for ancestry estimation.  We provide two sets of plots.  

-   GAP:  This visualization the proportion of each ancestry in individual samples.
-   LAP:  This visualization shows most probable posterior ancestry by regions of the chromosome.

### Module 9: PCA

We construct principal components to be used in future analysis (especially for correction with respect to ancestry).

### Module 10: population stratification

We separate the samples into individual Plink files based on their most probable posterior ancestries determined in Module 7: rfmix.

### Module 11: Subpopulation QC

We also provide individual QC steps stratified by population.

## Optional Pipeline Output

-   Gather information from log files created by standardized steps and create an automated report with tables and figures regarding each step.



## Docker ![Docker](https://img.shields.io/badge/docker-%230db7ed.svg?style=for-the-badge&logo=docker&logoColor=white)
Still in development

## Individual Module Usage

### Module 1: Crossmap

Requirements:

-	CrossMap.py
-	GRCh37_to_GRCh38.chain.gz

I recommend converting to Plink and performing these steps:

1. `plink --bfile file_stem --merge-x no-fail --make-bed --out prep1`
2. `plink --bfile prep1 --recode --output 'MT' --out prep2`
3. `awk '{print $1, $4-1, $4, $2}' prep2.map > prep.bed`
4. `python /path_to_crossmap_repo/CrossMap/CrossMap.py /path_to_crossmap_repo/CrossMap/GRCh37_to_GRCh38.chain.gz prep.bed study_lifted.bed`
5. `awk '{print $4}' study_lifted.bed > updated.snp`
6. `awk '{print $4, $3}' study_lifted.bed > updated.position`
7. `awk '{print $4, $1}' study_lifted.bed > updated.chr`
8. `plink --file prep2 --extract updated.snp --make-bed --out result1`
9. `plink --bfile result1 --update-map updated.position --make-bed --out result2`
10. `plink --bfile result2 --update-chr updated.chr --make-bed --out result3`
11. `plink --bfile result3 --recode --make-bed --out study_lifted`


### Module 2: GenotypeHarmonizer

Requirements:

-	GenotypeHarmonizer.jar
-	ALL.hgdp1kg.filtered.SNV_INDEL.38.phased.shapeit5

Make sure the data is in build GRCh38 first.  Subsequent steps:

1. `for chr in {1..22} X Y; do 
plink --bfile study_lifted --chr $chr --make-bed --out lifted.chr${chr};  done`
2. `java -Xmx8g -jar /path_to_genotype_harmonizer/GenotypeHarmonizer/GenotypeHarmonizer.jar --input lifted.chr${CHR} --inputType PLINK_BED --ref ALL.hgdp1kg.filtered.SNV_INDEL.38.phased.shapeit5 --refType VCF --keep --output lifted.chr${CHR}.aligned`

### Module 3: Initial QC

Initial QC recommends missing by individual and missing by genotype filters of 2 percent.

1. `plink --bfile file_stem --geno 0.02 --make-bed QC1`
2. `plink --bfile QC1 --mind 0.02 --make-bed --out QC2`

### Module 4: Relatedness

We need to check whether the FIDs are all 0.

`all_fids_zero=$(awk '{if ($1 != 0) exit 1}' temp.fam && echo "yes" || echo "no")`
`if [ "$all_fids_zero" == "yes" ]; then
    echo "All FIDs are zero. Using IID as FID for KING."
    # Save mapping to restore later
    awk '{print $2, $1}' temp.fam > original_fid_map.txt
    # Replace FID with IID
    awk '{print $2,$2,$3,$4,$5,$6}' temp.fam > kin.fam
else
    echo "FIDs are not all zero. Using original FID/IID."
    cp temp.fam kin.fam
fi`

We run kinship:

Requirements

-	Need to have king executable file
-	Recommend having our pca_ir_pipeline.R file readily available

`/path_to_repo/king -b QC_Initial.bed --kinship --prefix kinships`

To run PCAIR and PCRelate

`Rscript /path_to_repo/src/pca_ir_pipeline.R $WORK_DIRECTORY $FILE_NAME`

### Module 5: Standard QC

Here are the recommended steps for Standard QC (must be on unrelated individuals)

1. `plink --bfile study_unrelated --geno 0.1 --make-bed --out QC1`
2. `plink --bfile QC1 --mind 0.1 --make-bed --out QC2`
3. `plink --bfile QC2 --geno 0.02 --make-bed --out QC3`
4. `plink --bfile QC3 --mind 0.02 --make-bed --out QC4`
5. `plink --bfile QC4 --missing`
6. `plink --bfile QC4 --check-sex`
7. `grep 'PROBLEM' plink.sexcheck| awk '{print$1,$2}'>sex_discrepancy.txt`
8. `plink --bfile QC4 --remove sex_discrepancy.txt --make-bed --out QC5`
9. `plink --bfile QC5 --freq --out MAF_check`
10. `plink --bfile --maf 0.01 --make-bed --out QC6`
11. `plink --bfile QC6 --hardy`
12. `awk '{ if ($9 <0.00001) print $0 }' plink.hwe>plinkzoomhwe.hwe`
13. `plink --bfile QC6 --hwe 1e-6 --make-bed --out QC7a`
14. `plink --bfile QC7a --hwe 1e-10 --hwe-all --make-bed --out QC7`
15. `plink --bfile QC7 --exclude /path_to_github_repo/data/inversion.txt --range --indep-pairwise 50 5 0.2 --out indepSNP`
16. `plink --bfile QC7 --extract indepSNP.prune.in --het --out R_check`
17. `Rscript --no-save /path_to_github_repo/src/heterozygosity_outliers_list.R`
18. `sed 's/"// g' fail-het-qc.txt | awk '{print$1, $2}'> het_fail_ind.txt`
19. `plink --bfile QC7 --make-bed --out QC_final`

### Module 6: Phasing

Requirements
-	shapeit4.2
-	chr${CHR}.b38.gmap.gz

Prepare QC data set `QC8` and break it apart by chromosome and recode to vcf

1. `plink --bfile QC8 --chr $CHR --recode vcf --out study.chr${CHR}`
2. `bgzip -c study.chr${CHR}.vcf > study.chr${CHR}.vcf.gz`
3. `bcftools index -f study.chr${CHR}.vcf.gz`

Run phasing through shapeit4

4. `shapeit4.2 --input study.chr${CHR}.vcf.gz --map chr${CHR}.b38.gmap.gz --region ${CHR} --output study.chr${CHR}.phased.vcf.gz --thread`

### Module 7: rfmix

RFMIX requires the files be phased and in vcf.gz form.  
Requirements:
-	`hg38_phased.vcf.gz`
-	`super_population_map_file.txt`
-	`genetic_map_hg38.txt`

To run rfmix:
`rfmix -f study.chr${CHR}.phased.vcf.gz -r hg38_phased.vcf.gz -m super_population_map_file.txt -g genetic_map_hg38.txt -o ancestry_chr${CHR} --chromosome=$CHR`

### Module 8: ancestry plots

This module provides visualization for ancestry estimation.  We provide two sets of plots.  

-   GAP:  This visualization the proportion of each ancestry in individual samples.

Need to prepare files like this:

```
for chr in {1..22}; do
    input_file="$WORK/rfmix/ancestry_chr${chr}.rfmix.Q"
    for ind in $(seq 3 "$n_rfmix_rows"); do
        individual_index=$((ind - 2))
        output_file="$WORK/visualization/ancestry${individual_index}_chr${chr}.rfmix.Q"
        sed -n -e "1p" -e "2p" -e "${ind}p" "$input_file" > "$output_file"
        if [ "$chr" -eq 1 ]; then
            sample_index=$((individual_index - 1))
            echo -e "${individual_index}\t${sample_ids[$sample_index]}" >> "$mapping_file"
        fi
    done
done
```
Pipeline for Global Ancestry Plots

```
python ${REF}/RFMIX2-Pipeline-to-plot/GAP/Scripts/RFMIX2ToBed4GAP.py --prefix $WORK/visualization/ancestry --output $WORK/visualization
python ${REF}/RFMIX2-Pipeline-to-plot/GAP/Scripts/BedToGap.py --input ancestry.bed --out ancestry_GAP.bed
input_bed="ancestry_GAP.bed"
gap_script="${REF}/RFMIX2-Pipeline-to-plot/GAP/Scripts/GAP.py"
output_dir="./GAP_individual_plots"
mkdir -p "$output_dir"
head -n 2 "$input_bed" > header.tmp
tail -n +3 "$input_bed" | cut -f1 | sort -u > individual_ids.txt
batch_num=1
ids_batch=()
while read -r id; do
    ids_batch+=("$id")
    if [ "${#ids_batch[@]}" -eq 10 ]; then
        batch_file="$output_dir/batch_${batch_num}.bed"
        batch_pdf="$output_dir/batch_${batch_num}.pdf"
        cat header.tmp > "$batch_file"
        for id_in_batch in "${ids_batch[@]}"; do
            grep -w "$id_in_batch" "$input_bed" >> "$batch_file"
        done
        echo "Generating $batch_pdf..."
        python "$gap_script" --input "$batch_file" --output "$batch_pdf"
        rm "$batch_file"
        ids_batch=()
        ((batch_num++))
    fi
done < individual_ids.txt
if [ "${#ids_batch[@]}" -gt 0 ]; then
    batch_file="$output_dir/batch_${batch_num}.bed"
    batch_pdf="$output_dir/batch_${batch_num}.pdf"
    cat header.tmp > "$batch_file"
    for id_in_batch in "${ids_batch[@]}"; do
        grep -w "$id_in_batch" "$input_bed" >> "$batch_file"
    done
    echo "Generating $batch_pdf..."
    python "$gap_script" --input "$batch_file" --output "$batch_pdf"
    rm "$batch_file"
fi
```

-   LAP:  This visualization shows most probable posterior ancestry by regions of the chromosome.

```
n_subs=${#sample_ids[@]}
for chr in {1..22}; do
    input_file="$WORK/rfmix/ancestry_chr${chr}.msp.tsv"
    for ind in $(seq 1 "$n_subs"); do
        ind1=$(( (2 * ind - 1) + 6 ))
        ind2=$(( (2 * ind) + 6 ))
        output_file="$WORK/visualization/ancestry${ind}_chr${chr}.msp.tsv"
        echo "Processing: $output_file (Columns: 1-6, $ind1, $ind2)"
        cut -f1-6,$ind1,$ind2 "$input_file" > "$output_file"
    done
done
python ${REF}/RFMIX2-Pipeline-to-plot/LAP/Scripts/RFMIX2ToBed.py --prefix $WORK/visualization/ancestry --output $WORK/visualization
mkdir -p $WORK/LAP_plots
for ind in $(seq 1 "$n_subs"); do
    input_file_1="$WORK/visualization/ancestry${ind}_hap1.bed"
    input_file_2="$WORK/visualization/ancestry${ind}_hap2.bed"
    output_file="$WORK/LAP_plots/ancestry${ind}.bed"
    output_LAP="$WORK/LAP_plots/ancestry${ind}.pdf"
    echo "Processing: $output_file"
    python ${REF}/RFMIX2-Pipeline-to-plot/LAP/Scripts/BedToLAP.py --bed1 "$input_file_1" --bed2 "$input_file_2" --out "$output_file"
    python ${REF}/RFMIX2-Pipeline-to-plot/LAP/Scripts/LAP.py -I "$output_file" -O "$output_LAP"
done
```
### Module 9: PCA

1. `sh commvar.sh hg38 phased study_unrelated refpref stupref`
2. `awk '{$6=1; print}' refpref.fam > refpref_recode.fam`
3. `mv refpref_recode.fam refpref.fam`
4. `awk '{$6=2; print}' stupref.fam > stupref_recode.fam`
5. `mv stupref_recode.fam stupref.fam`
6. `plink --bfile refpref --write-snplist --out ref_snps`
7. `plink --bfile stupref --extract ref_snps.snplist --make-bed --out stupref_common`
8. `plink --bfile refpref --extract ref_snps.snplist --make-bed --out refpref_common`
9. `echo stupref_common > mergelist.txt`
10. `plink --bfile stupref_common --biallelic-only strict --make-bed --out stupref_common_bi_tmp`
11. `plink --bfile refpref_common --biallelic-only strict --make-bed --out refpref_common_bi_tmp`
12. `plink --bfile stupref_common_bi_tmp --freq --out freq_study`
13. `plink --bfile refpref_common_bi_tmp --freq --out freq_ref`
14. `awk 'NR > 1 { print $2, $3, $4 }' freq_study.frq > study_alleles.txt`
15. `awk 'NR > 1 { print $2, $3, $4 }' freq_ref.frq > ref_alleles.txt`
16. `sort study_alleles.txt > study_alleles.sorted.txt`
17. `sort ref_alleles.txt > ref_alleles.sorted.txt`
18. `join -1 1 -2 1 study_alleles.sorted.txt ref_alleles.sorted.txt > joined_alleles.txt`
19. `awk '($2 == $4 && $3 == $5) || ($2 == $5 && $3 == $4)' joined_alleles.txt | cut -d' ' -f1 > consistent_snps.txt`
20. `plink --bfile stupref_common_bi_tmp --extract consistent_snps.txt --make-bed --out stupref_common_bi`
21. `plink --bfile refpref_common_bi_tmp --extract consistent_snps.txt --make-bed --out refpref_common_bi`
22. `echo "refpref_common_bi" > merge_list.txt`
23. `plink --bfile stupref_common_bi --merge-list merge_list.txt --make-bed --out merged_common_bi --allow-no-sex`
24. `plink --bfile merged_common_bi --pca --out merged_dataset_pca --allow-no-sex`

### Module 10: population stratification

We separate the samples into individual Plink files based on their most probable posterior ancestries determined in Module 7: rfmix.
### Subpopulation QC

We can also provide individual QC steps stratified by population.


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
