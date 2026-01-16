#!/bin/bash

INPUT=$1
OUT=$2
REF=$3

source /projects/standard/gdc/public/envs/load_miniconda3-2.sh
conda activate gdcPipeline

mkdir -p $OUT

sh scripts/commvar.sh ${REF}/rfmix_ref/hg38_phased $INPUT $OUT/refpref $OUT/stupref

# For controls (reference panel)
awk '{$6=1; print}' $OUT/refpref.fam > $OUT/refpref_recode.fam
mv $OUT/refpref_recode.fam $OUT/refpref.fam

# For cases (study individuals)
awk '{$6=2; print}' $OUT/stupref.fam > $OUT/stupref_recode.fam
mv $OUT/stupref_recode.fam $OUT/stupref.fam

plink --bfile $OUT/refpref --write-snplist --out $OUT/ref_snps
plink --bfile $OUT/stupref --extract $OUT/ref_snps.snplist --make-bed --out $OUT/stupref_common
plink --bfile $OUT/refpref --extract $OUT/ref_snps.snplist --make-bed --out $OUT/refpref_common

echo $OUT/stupref_common > $OUT/mergelist.txt

plink --bfile $OUT/stupref_common --biallelic-only strict --make-bed --out $OUT/stupref_common_bi_tmp
plink --bfile $OUT/refpref_common --biallelic-only strict --make-bed --out $OUT/refpref_common_bi_tmp

# Step 1: Compute allele frequencies
plink --bfile $OUT/stupref_common_bi_tmp --freq --out $OUT/freq_study
plink --bfile $OUT/refpref_common_bi_tmp --freq --out $OUT/freq_ref

# Step 2: Extract variant ID and alleles
awk 'NR > 1 { print $2, $3, $4 }' $OUT/freq_study.frq > $OUT/study_alleles.txt
awk 'NR > 1 { print $2, $3, $4 }' $OUT/freq_ref.frq > $OUT/ref_alleles.txt

# Step 3: Sort and join
sort $OUT/study_alleles.txt > $OUT/study_alleles.sorted.txt
sort $OUT/ref_alleles.txt > $OUT/ref_alleles.sorted.txt
join -1 1 -2 1 $OUT/study_alleles.sorted.txt $OUT/ref_alleles.sorted.txt > $OUT/joined_alleles.txt

# Step 4: Keep SNPs with matching alleles
awk '($2 == $4 && $3 == $5) || ($2 == $5 && $3 == $4)' $OUT/joined_alleles.txt | cut -d' ' -f1 > $OUT/consistent_snps.txt

# Step 5: Filter and merge
plink --bfile $OUT/stupref_common_bi_tmp --extract $OUT/consistent_snps.txt --make-bed --out $OUT/stupref_common_bi
plink --bfile $OUT/refpref_common_bi_tmp --extract $OUT/consistent_snps.txt --make-bed --out $OUT/refpref_common_bi
echo "$OUT/refpref_common_bi" > $OUT/merge_list.txt
plink --bfile $OUT/stupref_common_bi --merge-list $OUT/merge_list.txt --make-bed --out $OUT/merged_common_bi --allow-no-sex

plink2 --bfile $OUT/merged_common_bi --pca approx --out $OUT/merged_dataset_pca --allow-no-sex --make-grm-bin
rm -f $OUT/*_tmp.* $OUT/snps_*.txt $OUT/intersect_snps.txt $OUT/merge_list.txt $OUT/ref_snps.*
