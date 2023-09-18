# Converting to VCF files so that LiftOver can be used
plink --bfile DATA  --recode vcf --out DATA

# Converting from VCF to BED (Browser Extensible Data)
grep -v '^#' DATA.vcf | awk -F '\t' '{print $1,$2,$2,$3}' > output.ucsc.bed

# Above worked 
# Now to add in the liftover part
# 
./liftOver output.ucsc.bed hg19ToHg38.over.chain.gz SMILES_GSA.conversions.ucsc.bed unMapped

# Assuming it works I'll need to convert it back into standard format...


