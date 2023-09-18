# Converting to VCF files so that CrossMap can be used

plink --bfile DATA  --recode vcf --out DATA 

# If fixed... it should be as easy as 
# Causing errors saying that it is unable to find CrossMap.py
# CrossMap.py vcf hg19ToHg38.over.chain.gz DATA.vcf hg38.fa out.hg38.vcf

# Converting back to bim/fam/bed format for the rest of the steps
plink --vcf out.hg38.vcf --make-bed --out hg38_DATA

# Assuming it works I'll need to convert it back into standard format...
