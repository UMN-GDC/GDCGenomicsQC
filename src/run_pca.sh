#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --time=48:00:00
#SBATCH -p agsmall
#SBATCH -o plinkPCA.out
#SBATCH -e plinkPCA.err
#SBATCH --job-name plinkPCA

WORK=$1
REF=$2
NAME=$3
path_to_repo=$4

module load plink

mkdir -p $WORK/PCA
cp $WORK/relatedness/study.$NAME.unrelated.* $WORK/PCA/
cd $WORK/PCA

sh ${path_to_repo}/src/commvar.sh ${REF}/rfmix_ref/hg38_phased study.$NAME.unrelated refpref stupref

# For controls (reference panel)
awk '{$6=1; print}' refpref.fam > refpref_recode.fam
mv refpref_recode.fam refpref.fam

# For cases (study individuals)
awk '{$6=2; print}' stupref.fam > stupref_recode.fam
mv stupref_recode.fam stupref.fam

plink --bfile refpref --write-snplist --out ref_snps
plink --bfile stupref --extract ref_snps.snplist --make-bed --out stupref_common
plink --bfile refpref --extract ref_snps.snplist --make-bed --out refpref_common

echo stupref_common > mergelist.txt

plink --bfile stupref_common --biallelic-only strict --make-bed --out stupref_common_bi_tmp
plink --bfile refpref_common --biallelic-only strict --make-bed --out refpref_common_bi_tmp

# Step 1: Compute allele frequencies
plink --bfile stupref_common_bi_tmp --freq --out freq_study
plink --bfile refpref_common_bi_tmp --freq --out freq_ref

# Step 2: Extract variant ID and alleles
awk 'NR > 1 { print $2, $3, $4 }' freq_study.frq > study_alleles.txt
awk 'NR > 1 { print $2, $3, $4 }' freq_ref.frq > ref_alleles.txt

# Step 3: Sort and join
sort study_alleles.txt > study_alleles.sorted.txt
sort ref_alleles.txt > ref_alleles.sorted.txt
join -1 1 -2 1 study_alleles.sorted.txt ref_alleles.sorted.txt > joined_alleles.txt

# Step 4: Keep SNPs with matching alleles
awk '($2 == $4 && $3 == $5) || ($2 == $5 && $3 == $4)' joined_alleles.txt | cut -d' ' -f1 > consistent_snps.txt

# Step 5: Filter and merge
plink --bfile stupref_common_bi_tmp --extract consistent_snps.txt --make-bed --out stupref_common_bi
plink --bfile refpref_common_bi_tmp --extract consistent_snps.txt --make-bed --out refpref_common_bi
echo "refpref_common_bi" > merge_list.txt
plink --bfile stupref_common_bi --merge-list merge_list.txt --make-bed --out merged_common_bi --allow-no-sex

plink --bfile merged_common_bi --pca --out merged_dataset_pca --allow-no-sex
rm -f *_tmp.* snps_*.txt intersect_snps.txt merge_list.txt ref_snps.*


# Step 6. UMAP reduction
# [ ] Need to make sure the filepaths make sense
#%% UMAP
Rscript umap.R --eigens merged_dataset_pca --out merged_dataset_umap \
  --npc 50 --neighbors 50 \
  --threads 1 \
  --ncoords 3 \
  --seed $RANDOM

# Step 7. Variational autoencoder
#recode plink bfile to vcf (prepping for popVAE)
plink2 --bfile merged_common_bi --recode vcf-iid --out merged_common_bi

# [ ] Still need to make sure the environment for this is compatible
/home/gdc/public/envs/popvae/bin/python /home/gdc/public/popvae/popvae.py \
    --infile merged_common_bi.vcf \
    --out merged_common_vae
