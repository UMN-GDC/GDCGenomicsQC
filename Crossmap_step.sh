# Adding in the virual environment where CrossMap is installed
source /home/faird/shared/code/external/envs/miniconda3/load_miniconda3.sh
conda activate GDC_pipeline

# Calling CrossMap
CrossMap.py vcf hg19ToHg38.over.chain.gz vcf_data hg38.fa out.hg38.vcf

conda deactivate # Get out of GDC_pipeline env.
conda deactivate # Get out of DCAN_base env.

