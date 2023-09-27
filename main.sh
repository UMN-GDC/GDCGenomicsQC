#!/bin/bash -l        
#SBATCH --time=4:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10g
#SBATCH --tmp=10g
#SBATCH -p agsmall
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=x500@umn.edu 

# Default values
INFILE="SMILES"
OUTDIR="path/to/output/directory"

# Check if directory exists, and create it if not
if [ ! -d "$OUTDIR" ]; then
    mkdir -p "$OUTDIR"
fi

# Function to display help message
show_help() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  --infile <prefix>       Specify the input file prefix (default: SMILES)"
    echo "  --outdir <directory>    Specify the output directory (default: /path/to/default/output/directory)"
    echo "  -h, --help              Show this help message"
    exit 1
}


# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --infile)
            INFILE="$2"
            shift 2
            ;;
        --outdir)
            OUTDIR="$2"
            shift 2
            ;;
        -h|--help)
            show_help
            ;;
        *)
            echo "Unknown argument: $1"
            show_help
            ;;
    esac
done

TEMPDIR=${OUTDIR}/temp/
if [ ! -d "$TEMPDIR" ]; then
    mkdir -p "$TEMPDIR"
fi



module load plink
module load python3/3.9.3_anaconda2021.11_mamba
module load R/4.2.2-openblas



#### Start of QC_steps_pt1.R conversion
echo "(step 1 of QC steps)"
# plink --bfile ${TEMPDIR}GSA --keep non_white_ids_b.txt --make-bed --out ${TEMPDIR}GSA_a

echo "(Step 1.5 of QC steps) removing samples with missingness more than 10%"
#This step is necessary only for the non_white abcd group 
# plink --bfile ${TEMPDIR}GSA_a --mind 0.1 --make-bed --out ${TEMPDIR}GSA_b
# Brought the sample size down to 1066, and 1402787 variants

echo "(step 2 of QC steps) remove markers where 10% or higher of samples have missing data"
#  at these markers 
# plink --bfile ${INFILE} --geno 0.1 --make-bed --out ${TEMPDIR}GSA_1

echo "(Step 3 of QC) Initial sample missing filtering"
plink --bfile ${INFILE} --mind 0.1 --make-bed --out ${TEMPDIR}GSA_2

echo "(Step 4 of QC) Ultimate marker missing filtering"
plink --bfile ${TEMPDIR}GSA_2 --geno 0.02 --make-bed --out ${TEMPDIR}GSA_3

echo "(Step 5 of QC) Ultimate sample missing filtering"
plink --bfile ${TEMPDIR}GSA_3 --mind 0.02 --make-bed --out ${TEMPDIR}GSA_4

echo "(Step 6 of QC) Check for sex discrepancies"
plink --bfile ${TEMPDIR}GSA_4 --check-sex
grep 'PROBLEM' plink.sexcheck| awk '{print$1,$2}'>sex_discrepancy.txt
plink --bfile ${TEMPDIR}GSA_4 --remove sex_discrepancy.txt --make-bed --out ${TEMPDIR}GSA_5

echo "(step 7 of QC) Minor allele frequency filtering"
plink --bfile ${TEMPDIR}GSA_4 --freq --out MAF_check # Allows generation of MAF distribution
plink --bfile ${TEMPDIR}GSA_4 --maf 0.01 --make-bed --out ${TEMPDIR}GSA_6

## Including below for use in src/QC_report.R ##
plink --bfile ${TEMPDIR}GSA_4 --hardy


echo "(step 8 of QC) HWE filter. Selecting SNPs with HWE p-value below 0.00001, allows to zoom in on strongly deviating SNPs."
awk '{ if ($9 <0.00001) print $0 }' plink.hwe>plinkzoomhwe.hwe
## ##

plink --bfile ${TEMPDIR}GSA_6 --hwe 1e-6 --make-bed --out ${TEMPDIR}GSA_7a
plink --bfile ${TEMPDIR}GSA_7a --hwe 1e-10 --hwe include-nonctrl --make-bed --out ${TEMPDIR}GSA_7

echo "(Step 8 of QC) Heterozygosity check"
plink --bfile ${TEMPDIR}GSA_7 --exclude data/inversion.txt --range --indep-pairwise 50 5 0.2 --out indepSNP
plink --bfile ${TEMPDIR}GSA_7 --extract indepSNP.prune.in --het --out R_check

# The following code generates a list of individuals who deviate more than 3 standard deviations from the heterozygosity rate mean.
Rscript --no-save src/heterozygosity_outliers_list.R 
sed 's/"// g' fail-het-qc.txt | awk '{print$1, $2}'> het_fail_ind.txt 
plink --bfile ${TEMPDIR}GSA_7 --remove het_fail_ind.txt --make-bed --out ${TEMPDIR}GSA_8

echo "(Step 9 of QC) Cryptic Relatedness check "
#Excluding individuals with parents in the study (only founders)

## Included to see if this fixes visualization issues ##
# Check for relationships between individuals with a pihat > 0.2.
plink --bfile ${TEMPDIR}GSA_8 --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2

# The following commands will visualize specifically these parent-offspring relations, using the z values. 
awk '{ if ($8 >0.9) print $0 }' pihat_min0.2.genome>zoom_pihat.genome
## ##

plink --bfile ${TEMPDIR}GSA_8 --filter-founders --make-bed --out ${TEMPDIR}GSA_9a

#Looking at individuals with pi_hat > 0.2
plink --bfile ${TEMPDIR}GSA_9a --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2_in_founders

#Remove individual of the pair that has pi_hat > 0.2, to generate the call rate list
plink --bfile ${TEMPDIR}GSA_9a --missing

# Generates PDF QC_report 
Rscript --no-save src/QC_report.R 

Rscript src/lower_pihat_list_generator.R
plink --bfile ${TEMPDIR}GSA_9a --remove 0.2_low_call_rate_pihat.txt --make-bed --out ${TEMPDIR}GSA_10

echo "QC steps are done!"

