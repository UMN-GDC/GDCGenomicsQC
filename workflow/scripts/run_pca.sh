#!/bin/bash

INPUT=$1
STAGE=$2
REF=$3

mkdir -p $STAGE/intermediates

sh scripts/commvar.sh ${REF}/rfmix_ref/hg38_phased $INPUT $STAGE/intermediates/refpref $STAGE/intermediates/stupref

# For controls (reference panel)
awk '{$6=1; print}' $STAGE/intermediates/refpref.fam > $STAGE/intermediates/refpref_recode.fam
mv $STAGE/intermediates/refpref_recode.fam $STAGE/intermediates/refpref.fam

# For cases (study individuals)
awk '{$6=2; print}' $STAGE/intermediates/stupref.fam > $STAGE/intermediates/stupref_recode.fam
mv $STAGE/intermediates/stupref_recode.fam $STAGE/intermediates/stupref.fam

plink --bfile $STAGE/intermediates/refpref --write-snplist --out $STAGE/ref_snps
plink --bfile $STAGE/intermediates/stupref --extract $STAGE/ref_snps.snplist --make-bed --out $STAGE/intermediates/stupref_common
plink --bfile $STAGE/intermediates/refpref --extract $STAGE/ref_snps.snplist --make-bed --out $STAGE/intermediates/refpref_common

echo $STAGE/intermediates/stupref_common > $STAGE/mergelist.txt

plink --bfile $STAGE/intermediates/stupref_common --biallelic-only strict --make-bed --out $STAGE/intermediates/stupref_common_bi_tmp
plink --bfile $STAGE/intermediates/refpref_common --biallelic-only strict --make-bed --out $STAGE/intermediates/refpref_common_bi_tmp

# Step 1: Compute allele frequencies
plink --bfile $STAGE/intermediates/stupref_common_bi_tmp --freq --out $STAGE/freq_study
plink --bfile $STAGE/intermediates/refpref_common_bi_tmp --freq --out $STAGE/freq_ref

# Step 2: Extract variant ID and alleles
awk 'NR > 1 { print $2, $3, $4 }' $STAGE/freq_study.frq > $STAGE/study_alleles.txt
awk 'NR > 1 { print $2, $3, $4 }' $STAGE/freq_ref.frq > $STAGE/ref_alleles.txt

# Step 3: Sort and join
sort $STAGE/study_alleles.txt > $STAGE/study_alleles.sorted.txt
sort $STAGE/ref_alleles.txt > $STAGE/ref_alleles.sorted.txt
join -1 1 -2 1 $STAGE/study_alleles.sorted.txt $STAGE/ref_alleles.sorted.txt > $STAGE/joined_alleles.txt

# Step 4: Keep SNPs with matching alleles
awk '($2 == $4 && $3 == $5) || ($2 == $5 && $3 == $4)' $STAGE/joined_alleles.txt | cut -d' ' -f1 > $STAGE/consistent_snps.txt

# Step 5: Filter and merge
plink --bfile $STAGE/intermediates/stupref_common_bi_tmp --extract $STAGE/consistent_snps.txt --make-bed --out $STAGE/intermediates/stupref_common_bi
plink --bfile $STAGE/intermediates/refpref_common_bi_tmp --extract $STAGE/consistent_snps.txt --make-bed --out $STAGE/intermediates/refpref_common_bi
echo "$STAGE/intermediates/refpref_common_bi" > $STAGE/merge_list.txt
plink --bfile $STAGE/intermediates/stupref_common_bi --merge-list $STAGE/merge_list.txt --make-bed --out $STAGE/merged_common_bi --allow-no-sex

# Compute PCs
plink2 --bfile $STAGE/merged_common_bi \
       --freq counts \
       --pca approx allele-wts vcols=chrom,ref,alt \
       --out $STAGE/merged_dataset_pca \
       --make-grm-bin \
       --allow-no-sex

# then project related onto those PCs if there are unrelated
if [ -f "${STAGE}/../02-relatedness/related.bed" ]; then
  echo "Project unrelated sample onto the PCs."
  plink2 --bfile ${STAGE}/../02-relatedness/related \
         --read-freq $STAGE/merged_dataset_pca.acount \
         --score $STAGE/merged_dataset_pca.eigenvec.allele 2 5 header-read no-mean-imputation \
                 variance-standardize \
         --score-col-nums 6-15 \
         --out $STAGE/merged_dataset_pca_related
fi
