#!/bin/bash

show_help() {
  echo "Usage: $0 [options]"
  echo "Options:"
  echo "  --FILE <file_name>  Specify the file to process. Must be a .bed file"
  echo "  --help              Display this help message."
}

# Check for command line arguments
if [ $# -eq 0 ]; then
  echo "No arguments provided. Use --help for usage information."
  exit 1
fi

# Loop through the command line arguments
while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
    --FILE)
      FILE="$2"
      shift 2 # Consume both the flag and its value
      ;;
    --help)
      show_help
      exit 0
      ;;
    *)
      echo "Unrecognized option: $key"
      show_help
      exit 1
      ;;
  esac
done

# Just in case they included the file extension
FILE=${FILE%*}


echo "(step 1 of QC steps) Matching data to NIH's GRCh38 genome build"

# Lift-over (if data's genome build is GRCh37 or earlier): adding in CrossMap 
# Converting to VCF for process
plink --bfile $FILE --recode vcf --out vcf_data
./Crossmap_step.sh
plink --vcf out.hg38.vcf --make-bed --out hg38_${FILE}

# Update strand orientation and flip alleles as needed:
# Note: SMILES-GSA data was genotyped using the chip array GSA-24v2-0
wget https://www.well.ox.ac.uk/~wrayner/strand/GSA-24v2-0_A2-b38-strand.zip
wget https://www.well.ox.ac.uk/~wrayner/strand/update_build.sh
unzip GSA-24v2-0_A2-b38-strand.zip && rm GSA-24v2-0_A2-b38-strand.zip
./update_build.sh hg38_${FILE} GSA-24v2-0_A2-b38.strand ${FILE}_hg38NIH

echo "(step 2 of QC steps) remove markers where 10% or higher of samples have missing data"
#  at these markers 
#plink --bfile ${FILE}_hg38NIH --geno 0.1 --make-bed --out ${FILE}_1

echo "(Step 3 of QC) Initial sample missing filtering"
plink --bfile ${FILE}_1 --mind 0.1 --make-bed --out ${FILE}_2

echo "(Step 4 of QC) Ultimate marker missing filtering"
plink --bfile ${FILE}_2 --geno 0.02 --make-bed --out ${FILE}_3

echo "(Step 5 of QC) Ultimate sample missing filtering"
plink --bfile ${FILE}_3 --mind 0.02 --make-bed --out ${FILE}_4

echo "(Step 6 of QC) Check for sex discrepancies"
plink --bfile ${FILE}_4 --check-sex
#To be added: we should modify the below grep code to include only non-ambiguous reported gender or ambiguous genotype-gender
#For missing reported gender, prompt if user want to impute using genotype-gender
grep 'PROBLEM' plink.sexcheck| awk '{print$1,$2}'>sex_discrepancy.txt

plink --bfile ${FILE}_4 --remove sex_discrepancy.txt --make-bed --out ${FILE}_5

echo "(step 7 of QC) Minor allele frequency filtering"
plink --bfile ${FILE}_4 --freq --out MAF_check # Allows generation of MAF distribution
plink --bfile ${FILE}_4 --maf 0.01 --make-bed --out ${FILE}_6

## Including below for use in src/QC_report.R ##
plink --bfile ${FILE}_4 --hardy


echo "(step 8 of QC) HWE filter. Selecting SNPs with HWE p-value below 0.00001, allows to zoom in on strongly deviating SNPs."
awk '{ if ($9 <0.00001) print $0 }' plink.hwe>plinkzoomhwe.hwe
## ##

plink --bfile ${FILE}_6 --hwe 1e-6 --make-bed --out ${FILE}_7a
plink --bfile ${FILE}_7a --hwe 1e-10 --hwe-all --make-bed --out ${FILE}_7

echo "(Step 8 of QC) Heterozygosity check"
plink --bfile ${FILE}_7 --exclude data/inversion.txt --range --indep-pairwise 50 5 0.2 --out indepSNP
plink --bfile ${FILE}_7 --extract indepSNP.prune.in --het --out R_check

# The following code generates a list of individuals who deviate more than 3 standard deviations from the heterozygosity rate mean.
Rscript --no-save src/heterozygosity_outliers_list.R 
sed 's/"// g' fail-het-qc.txt | awk '{print$1, $2}'> het_fail_ind.txt 
plink --bfile ${FILE}_7 --remove het_fail_ind.txt --make-bed --out ${FILE}_8

echo "(Step 9 of QC) Cryptic Relatedness check "
#Excluding individuals with parents in the study (only founders)

## Included to see if this fixes visualization issues ##
# Check for relationships between individuals with a pihat > 0.2.
plink --bfile ${FILE}_8 --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2

# The following commands will visualize specifically these parent-offspring relations, using the z values. 
awk '{ if ($8 >0.9) print $0 }' pihat_min0.2.genome>zoom_pihat.genome
## ##

plink --bfile ${FILE}_8 --filter-founders --make-bed --out ${FILE}_9a

#Looking at individuals with pi_hat > 0.2
plink --bfile ${FILE}_9a --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2_in_founders

#Remove individual of the pair that has pi_hat > 0.2, to generate the call rate list
plink --bfile ${FILE}_9a --missing

#To be added: we need a script to remove only 1 individual of the pair that has pi_hat > 0.2.
#This is like a covariance matrix where we have to remove one of the pair that correlate more than 0.2
#The problem being the order of removal since 1 sample might be related to several other samples, generating several lines
#in the file pihat_min0.2_in_founders

Rscript src/lower_pihat_list_generator.R
plink --bfile ${FILE}_9a --remove 0.2_low_call_rate_pihat.txt --make-bed --out ${FILE}_10

echo "(Step 10 of QC) Principle Component Analysis"
#We use fraposa to perform pca on the current data and a reference data 1000G with known population
#We superimpose the pcs of the current data onto the reference data to compare and predict the current data's population
git clone https://github.com/daviddaiweizhang/fraposa.git
mv ./fraposa/*.* ./ && rm -R fraposa
./commvar.sh 1000G ${FILE}_10 1000G_comm ${FILE}_11
./fraposa_runner.py --stu_filepref ${FILE}_11 1000G_comm
./predstupopu.py 1000G_comm ${FILE}_11
./plotpcs.py 1000G_comm ${FILE}_11

# Generates PDF QC_report 
Rscript --no-save src/QC_report.R 

echo "QC steps are done!"
