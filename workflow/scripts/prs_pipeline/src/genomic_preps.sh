#!/bin/bash

# --- 1. CONFIGURATION ---
# Paths to your container and data
CONTAINER="/home/ood-baron063/prs_pipeline/prsv2_latest.sif"
RAW_DATA="/shared/release/abcd/abcd/concatenated/genetics/genotype_microarray/smokescreen/merged_chroms"
PHENO_FILE="/home/ood-baron063/R25_testing/ready_pheno.txt"                  # FID IID Value
GENDER_FILE="/home/ood-baron063/R25_testing/ready_sex.txt"                # FID IID GenderCode
UPDATED_PREFIX="abcd_updated_all"       # Output name after updating pheno/sex

awk '{$1=$1; print}' OFS='\t' "$GENDER_FILE" > gender_tab_fixed.txt
awk '{$1=$1; print}' OFS='\t' "$PHENO_FILE" > pheno_tab_fixed.txt


# --- 2. DEFINE THE RUN COMMAND ---
# This helper alias saves typing the container prefix every time
PLINK_EXEC="apptainer exec $CONTAINER plink"
PLINK2_EXEC="apptainer exec $CONTAINER plink2"

echo "Starting PLINK Pipeline..."

# --- STEP 1: UPDATE PHENOTYPE AND GENDER ---
echo "Step 1: Updating phenotype and gender information..."
$PLINK_EXEC --bfile $RAW_DATA \
            --update-sex gender_tab_fixed.txt \
            --pheno pheno_tab_fixed.txt \
            --make-bed \
            --out $UPDATED_PREFIX


# --- STEP 2: SPLIT BY ANCESTRY AND QUALITY CONTROL (PLINK 1.9) ---
echo "Step 2a: Extracting and filtering EUR ancestry..."
$PLINK_EXEC --bfile $UPDATED_PREFIX \
            --keep /home/ood-baron063/R25_testing/keep_eur.txt \
            --geno 0.05 \
            --mind 0.05 \
            --make-bed \
            --out abcd_EUR_filtered

echo "Step 2b: Extracting and filtering AFR ancestry..."
$PLINK_EXEC --bfile $UPDATED_PREFIX \
            --keep /home/ood-baron063/R25_testing/keep_afr.txt \
            --geno 0.05 \
            --mind 0.05 \
            --make-bed \
            --out abcd_AFR_filtered


# --- STEP 3: CALCULATE ALLELE FREQUENCIES (PLINK 2) ---
echo "Step 3a: Calculating allele frequencies for EUR..."
$PLINK2_EXEC --bfile abcd_EUR_filtered \
             --freq \
             --make-bed \
             --out abcd_EUR_final

echo "Step 3b: Calculating allele frequencies for AFR..."
$PLINK2_EXEC --bfile abcd_AFR_filtered \
             --freq \
             --make-bed \
             --out abcd_AFR_final

echo "All processes ran."