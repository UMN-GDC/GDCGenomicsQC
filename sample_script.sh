#!/bin/bash -l        
#SBATCH --time=4:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10g
#SBATCH --tmp=10g
#SBATCH -p agsmall
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=baron063@umn.edu 

module load R
module load plink

# Conversion to vcf format for CrossMap
Crossmap_step.sh
#Running CrossMap through Python
### Put script here ###

#(step 1 of QC steps)
# plink --bfile SMILES_GSA --keep non_white_ids_b.txt --make-bed --out SMILES_GSA_a

#(Step 1.5 of QC steps) removing samples with missingness more than 10%
#This step is necessary only for the non_white abcd group 
# plink --bfile SMILES_GSA_a --mind 0.1 --make-bed --out SMILES_GSA_b


#(step 2 of QC steps) remove markers where 10% or higher of samples have missing data 
#  at these markers 
# plink --bfile SMILES_GSA --geno 0.1 --make-bed --out SMILES_GSA_1

#(Step 3 of QC) Initial sample missing filtering
plink --bfile SMILES_GSA --mind 0.1 --make-bed --out SMILES_GSA_2

#(Step 4 of QC) Ultimate marker missing filtering
plink --bfile SMILES_GSA_2 --geno 0.02 --make-bed --out SMILES_GSA_3

#(Step 5 of QC) Ultimate sample missing filtering
plink --bfile SMILES_GSA_3 --mind 0.02 --make-bed --out SMILES_GSA_4

# #(Step 6 of QC) Check for sex discrepancies
plink --bfile SMILES_GSA_4 --check-sex
grep 'PROBLEM' plink.sexcheck| awk '{print$1,$2}'>sex_discrepancy.txt
plink --bfile SMILES_GSA_4 --remove sex_discrepancy.txt --make-bed --out SMILES_GSA_5

#(step 7 of QC) Minor allele frequency filtering
plink --bfile SMILES_GSA_4 --freq --out MAF_check # Allows generation of MAF distribution
plink --bfile SMILES_GSA_4 --maf 0.01 --make-bed --out SMILES_GSA_6

## Including below for use in QC_report.R ##
plink --bfile SMILES_GSA_4 --hardy
# Selecting SNPs with HWE p-value below 0.00001, allows to zoom in on strongly deviating SNPs. 
awk '{ if ($9 <0.00001) print $0 }' plink.hwe>plinkzoomhwe.hwe
## ##

plink --bfile SMILES_GSA_6 --hwe 1e-6 --make-bed --out SMILES_GSA_7a
plink --bfile SMILES_GSA_7a --hwe 1e-10 --hwe-all --make-bed --out SMILES_GSA_7

#(Step 8 of QC) Heterozygosity check
plink --bfile SMILES_GSA_7 --exclude inversion.txt --range --indep-pairwise 50 5 0.2 --out indepSNP
plink --bfile SMILES_GSA_7 --extract indepSNP.prune.in --het --out R_check

# The following code generates a list of individuals who deviate more than 3 standard deviations from the heterozygosity rate mean.
Rscript --no-save heterozygosity_outliers_list.R 
sed 's/"// g' fail-het-qc.txt | awk '{print$1, $2}'> het_fail_ind.txt 
plink --bfile SMILES_GSA_7 --remove het_fail_ind.txt --make-bed --out SMILES_GSA_8

#(Step 9 of QC) Cryptic Relatedness check 
#Excluding individuals with parents in the study (only founders)

## Included to see if this fixes visualization issues ##
# Check for relationships between individuals with a pihat > 0.2.
plink --bfile SMILES_GSA_8 --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2

# The following commands will visualize specifically these parent-offspring relations, using the z values. 
awk '{ if ($8 >0.9) print $0 }' pihat_min0.2.genome>zoom_pihat.genome
## ##

plink --bfile SMILES_GSA_8 --filter-founders --make-bed --out SMILES_GSA_9a

#Looking at individuals with pi_hat > 0.2
plink --bfile SMILES_GSA_9a --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2_in_founders

#Remove individual of the pair that has pi_hat > 0.2, to generate the call rate list
plink --bfile SMILES_GSA_9a --missing

# Generates PDF QC_report 
Rscript --no-save QC_report.R 

Rscript lower_pihat_list_generator.R
plink --bfile SMILES_GSA_9a --remove 0.2_low_call_rate_pihat.txt --make-bed --out SMILES_GSA_10

#### QC steps are done!

