#!/bin/bash
# ----------------------------------------------------------
# Remove problematic SNPs from all per-chromosome PLINK files
# study.hbcd_rsid.lifted.chr${CHR}.aligned
# using exclusion list: study.hbcd_rsid.lifted.aligned1-merge.missnp
# ----------------------------------------------------------

MISSNP="study.hbcd_rsid.lifted.aligned1-merge.missnp"

# Ensure the missnp file exists
if [ ! -f "$MISSNP" ]; then
  echo "ERROR: Cannot find $MISSNP"
  exit 1
fi

# Loop over 1–22 + X + Y
for CHR in {1..22} X Y; do
  INFILE="study.hbcd_rsid.lifted.chr${CHR}.aligned"
  OUTFILE="study.hbcd_rsid.lifted.chr${CHR}.aligned.clean"

  if [ -f "${INFILE}.bed" ]; then
    echo "Processing chromosome ${CHR} ..."
    plink \
      --bfile "$INFILE" \
      --exclude "$MISSNP" \
      --make-bed \
      --out "$OUTFILE"
  else
    echo "WARNING: ${INFILE}.bed not found, skipping chromosome ${CHR}."
  fi
done

echo "All chromosomes processed. Cleaned files written as *.aligned.clean.{bed,bim,fam}"
