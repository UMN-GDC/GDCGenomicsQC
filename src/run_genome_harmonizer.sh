#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --time=12:00:00
#SBATCH -p agsmall
#SBATCH -o GenotypeHarmonizer.out
#SBATCH -e GenotypeHarmonizer.err
#SBATCH --job-name GenotypeHarmonizer

WORK=$1
REF=$2
NAME=$3
path_to_repo=$4
file_to_use=$5

echo "WORK: $WORK"
echo "REF: $REF"
echo "NAME: $NAME"
echo "path_to_repo: $path_to_repo"
echo "file_to_use: $file_to_use"

mkdir -p "$WORK/lifted"

# Split by chr (1–22, X, Y) after liftover
for chr in {1..22} X Y; do 
  plink \
    --bfile "${file_to_use}" \
    --chr "${chr}" \
    --make-bed \
    --out "$WORK/lifted/study.${NAME}.lifted.chr${chr}"
done

## Removed this from the above plink command to align with reference genome # --output-chr chrMT  

echo "Deleting extra files"
rm -f prep1.* prep2.* result1.* result2.* result3.* prep.bed updated.snp updated.position updated.chr

# Using genome harmonizer, update strand orientation and flip alleles according to the reference dataset.
echo "Begin autosomal harmonization"
mkdir -p "$WORK/aligned"

sbatch \
  --time=24:00:00 \
  --mem=15GB \
  --array=1-22 \
  --wait \
  -N1 \
  "${path_to_repo}/src/harmonizer_individual.job" \
  "${WORK}" \
  "${NAME}" \
  "${REF}"

mkdir -p "${WORK}/logs" "${WORK}/logs/errors" "${WORK}/logs/out"
mv -f "${WORK}"/*.out "${WORK}/logs/out/" 2>/dev/null || true
mv -f "${WORK}"/*.err "${WORK}/logs/errors/" 2>/dev/null || true

# Currently reference dataset does not have chrY for alignment, and ChrX has no match with study data
# Hence, we bring the unaligned ChrX and ChrY to the result folder, i.e. skipping alignment
cp "$WORK/lifted/study.${NAME}.lifted.chrX.bed" "$WORK/aligned/study.${NAME}.lifted.chrX.aligned.bed"
cp "$WORK/lifted/study.${NAME}.lifted.chrX.bim" "$WORK/aligned/study.${NAME}.lifted.chrX.aligned.bim"
cp "$WORK/lifted/study.${NAME}.lifted.chrX.fam" "$WORK/aligned/study.${NAME}.lifted.chrX.aligned.fam"

cp "$WORK/lifted/study.${NAME}.lifted.chrY.bed" "$WORK/aligned/study.${NAME}.lifted.chrY.aligned.bed"
cp "$WORK/lifted/study.${NAME}.lifted.chrY.bim" "$WORK/aligned/study.${NAME}.lifted.chrY.aligned.bim"
cp "$WORK/lifted/study.${NAME}.lifted.chrY.fam" "$WORK/aligned/study.${NAME}.lifted.chrY.aligned.fam"

# Parse GenotypeHarmonizer logs
"${path_to_repo}/src/genotype_harmonizer_log_reader.sh" "${WORK}/aligned"
## Creates genome_harmonizer_full_log.txt inside of the aligned directory

##############################################
#  MERGE BLOCK WITH 3+ ALLELE REMOVAL
##############################################

cd "${WORK}/aligned"

# Get all per-chromosome aligned BIMs (includes 1-22, X, Y)
temp_1=$(ls "study.${NAME}.lifted.chr"*.aligned.bim)
# Strip ".bim" to get the PLINK prefix list
array=(${temp_1//.bim/})

# First: attempt a merge to force PLINK to produce the -merge.missnp file
echo "Initial merge to detect multi-allelic variants..."
printf "%s\n" "${array[@]:1}" > mergelist_initial.txt

out_prefix="study.${NAME}.lifted.aligned1"

# Allow failure here so the script continues even if PLINK errors out
plink \
  --bfile "${array[0]}" \
  --merge-list mergelist_initial.txt \
  --allow-no-sex \
  --make-bed \
  --out "${out_prefix}" || true

missnp_file="${out_prefix}-merge.missnp"

if [ -f "${missnp_file}" ]; then
  echo "Multi-allelic / problematic variants detected in ${missnp_file}."
  echo "Excluding these variants from each chromosome (1–22, X, Y)..."

  # Loop over each per-chromosome dataset and drop those variants
  for chr_prefix in "${array[@]}"; do
    echo "  Cleaning ${chr_prefix} ..."
    plink \
      --bfile "${chr_prefix}" \
      --exclude "${missnp_file}" \
      --make-bed \
      --out "${chr_prefix}.noMulti"
  done

  # Rebuild the list of per-chromosome prefixes, now using the cleaned .noMulti sets
  temp_2=$(ls "study.${NAME}.lifted.chr"*.aligned.noMulti.bim)
  array_clean=(${temp_2//.bim/})

  # Now do the final merge with cleaned data
  echo "Final merge after excluding multi-allelic variants..."
  printf "%s\n" "${array_clean[@]:1}" > mergelist.txt

  plink \
    --bfile "${array_clean[0]}" \
    --merge-list mergelist.txt \
    --allow-no-sex \
    --make-bed \
    --out "${out_prefix}"

else
  echo "No -merge.missnp file found; proceeding with standard merge."
  # No multi-allelic issue detected → do a straight merge with original array
  printf "%s\n" "${array[@]:1}" > mergelist.txt

  plink \
    --bfile "${array[0]}" \
    --merge-list mergelist.txt \
    --allow-no-sex \
    --make-bed \
    --out "${out_prefix}"
fi

# Now proceed with your original two lines
plink \
  --bfile "${out_prefix}" \
  --split-x 'hg38' no-fail \
  --make-bed \
  --out "study.${NAME}.lifted.aligned"

echo "Done. Final merged + split-X dataset: study.${NAME}.lifted.aligned.{bed,bim,fam}"
