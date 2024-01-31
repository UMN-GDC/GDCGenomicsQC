#!/bin/bash

show_help() {
  echo "Usage: $0 [options]"
  echo "Options:"
  echo "  --FILE <file_name>    Specify the file to process. Must be a .bed file"
  echo "  --PATHTODATA          Specify the full path to where this file is stored"
  echo "  --PATHTOSTOREOUTPUTS  Specify the full path to where you would like the outputs of this pipeline to go"
  echo "  --OLDBLD <Y/N>        Specify if data's genome build is GRCh37 or earlier"
  echo "  --UPDTSO <Y/N>        Specify if you would like to update the strand orientation and flip alleles as necessary"
  echo "  --help                Display this help message."
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
        echo "File chosen is $FILE"
        shift 2 # Consume both the flag and its value
      ;;
    --PATHTODATA)
        path_to_data="$2"
        echo "Path provided is $path_to_data"
        shift 2
      ;;
    --PATHTOSTOREOUTPUTS)
        path_to_store_outputs="$2"
        echo "Path provided is $path_to_data"
        shift 2
      ;;   
    --help)
        show_help
        exit 0
      ;;
    --OLDBLD)
        OLDBLD="$2"
        case $OLDBLD in
          [yY]|[yY][eE][sS])
            echo "You have chosen to use Crossmap."
            OLDBLD="YES"
          ;;
          [nN]|[nN][oO])
            echo "You have chosen to not use Crossmap."
            OLDBLD="NO"
          ;;
          *) 
          echo "Please enter a valid option (Yes/No) when choosing to include the OLDBLD flag."
          exit 1
          ;;
        esac
        shift 2
      ;;
    --UPDTSO)
        UPDTSO="$2"
#        echo "You have entered $UPDTSO"
        case $UPDTSO in
          [yY]|[yY][eE][sS])
            echo "You have chosen to update the strand orientation and flip alleles as necessary."
            UPDTSO="YES"
          ;;
          [nN]|[nN][oO])
            echo "You have chosen to not update the strand orientation. "
            UPDTSO="NO"
          ;;
          *) 
            echo "Please enter a valid option (Yes/No) when choosing to include the UPDTSO flag."
            exit 1
          ;;
        esac
        shift 2
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

# Temporarily copying over the desired files into current directory
cp ${path_to_data}/${FILE}.fam .
cp ${path_to_data}/${FILE}.bim .
cp ${path_to_data}/${FILE}.bed .

# Making a directory for the log files, & this temporary step
mkdir -p ./tmp
mkdir -p ${path_to_store_outputs}/logs
mkdir -p ${path_to_store_outputs}/data

# For if they want to update the genome build to GRCh38
case "$OLDBLD" in
  [Y][E][S]) 
    echo "Matching data to NIH's GRCh38 genome build"
    FILE_B="cross"
    # Lift-over (if data's genome build is GRCh37 or earlier): adding in CrossMap 
#    plink --bfile $FILE --recode vcf --out vcf_data
#    ./src/Crossmap_step.sh
#    echo "${FILE_B}_${FILE}"
#    plink --vcf out.hg38.vcf --make-bed --out $FILE_B_${FILE}
#    exit 0
  ;;
  *) 
    FILE_B="same"
    plink --bfile ${FILE} --make-bed --out ${FILE_B}_${FILE}
#    echo "${FILE_B}_${FILE}"
#    exit 0
  ;;
esac


# Update strand orientation and flip alleles as needed:
# Note: SMILES-GSA data was genotyped using the chip array GSA-24v2-0
case "$UPDTSO" in
  [Y][E][S]) 
    echo "Updating strand orientation and flipping alleles as necessary"
    FILE_C="strand"
    wget https://www.well.ox.ac.uk/~wrayner/strand/GSA-24v2-0_A2-b38-strand.zip
    wget https://www.well.ox.ac.uk/~wrayner/strand/update_build.sh
    unzip GSA-24v2-0_A2-b38-strand.zip && rm GSA-24v2-0_A2-b38-strand.zip
    ./update_build.sh ${FILE_B}_${FILE} GSA-24v2-0_A2-b38.strand ${FILE_B}_${FILE_C}_${FILE} #This line gives the error "Permission denied"
#    echo "${FILE_B}_${FILE_C}_${FILE}"
#    exit 0
  ;;
  *) 
    FILE_C="same"
    plink --bfile ${FILE_B}_${FILE} --make-bed --out ${FILE_B}_${FILE_C}_${FILE}
#    echo "${FILE_B}_${FILE_C}_${FILE}"
#    exit 0
  ;;
esac


# Standard QC steps

echo "(step 2 of QC steps) remove markers where 10% or higher of samples have missing data"
#  at these markers 
plink --bfile ${FILE_B}_${FILE_C}_${FILE} --geno 0.1 --make-bed --out ${FILE}_1

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

Rscript src/lower_pihat_list_generator_v2.R
plink --bfile ${FILE}_9a --remove 0.2_low_call_rate_pihat.txt --make-bed --out ${FILE}_10

# Moving files that will be used by next steps
mv *.log ${path_to_store_outputs}/logs/
mv ${FILE}_10* ./tmp/

# Removing all intermediary steps
rm ${FILE}* 

# Putting main output back in this location
mv ./tmp/${FILE}_10* .

work_dir=$(pwd)
echo ${work_dir}

cp ${FILE}_10* ${path_to_store_outputs}/data/

echo "QC steps are done!"



