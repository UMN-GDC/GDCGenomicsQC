#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <sample_n> <cache_root> [seed]" >&2
  exit 1
fi

SAMPLE_N="$1"
CACHE_ROOT="$2"
SEED="${3:-42}"

PLINK2="${PLINK2:-/projects/standard/gdc/public/plink2}"
NPC="${NPC:-10}"
THREADS="${THREADS:-4}"
MEM_MB="${MEM_MB:-24000}"

declare -A INPUT_PREFIXES
INPUT_PREFIXES[AFR]="${AFR_INPUT_PREFIX:-/scratch.global/saonli/GDCGenomicsQC/CTSLEB/AFR/CTSLEB_AFR}"
INPUT_PREFIXES[EUR]="${EUR_INPUT_PREFIX:-/scratch.global/saonli/GDCGenomicsQC/CTSLEB/EUR/CTSLEB_EUR}"

mkdir -p "$CACHE_ROOT/N_${SAMPLE_N}"

for ANC in AFR EUR; do
  INPUT_PREFIX="${INPUT_PREFIXES[$ANC]}"
  PREFIX="$CACHE_ROOT/N_${SAMPLE_N}/${ANC}_spline"
  LOG="$CACHE_ROOT/N_${SAMPLE_N}/${ANC}_precompute.log"

  echo "===== $ANC / N=$SAMPLE_N ====="

  if [[ -s "${PREFIX}.grm.bin" && -s "${PREFIX}.grm.N.bin" && -s "${PREFIX}.grm.id" && -s "${PREFIX}.eigenvec" && -s "${PREFIX}_adjhe.eigenvec" ]]; then
    echo "Using existing cached files at ${PREFIX}" | tee "$LOG"
    continue
  fi

  mkdir -p "$(dirname "$PREFIX")"

  if [[ -s "${INPUT_PREFIX}.pgen" ]]; then
    "$PLINK2" --pfile "$INPUT_PREFIX" \
      --thin-indiv-count "$SAMPLE_N" --seed "$SEED" \
      --make-bed \
      --out "$PREFIX" \
      --threads "$THREADS" \
      --memory "$MEM_MB" > "$LOG" 2>&1
  elif [[ -s "${INPUT_PREFIX}.bed" ]]; then
    "$PLINK2" --bfile "$INPUT_PREFIX" \
      --thin-indiv-count "$SAMPLE_N" --seed "$SEED" \
      --make-bed \
      --out "$PREFIX" \
      --threads "$THREADS" \
      --memory "$MEM_MB" > "$LOG" 2>&1
  else
    echo "No .pgen or .bed input found for ${INPUT_PREFIX}" | tee -a "$LOG" >&2
    exit 1
  fi

  awk 'BEGIN {OFS="\t"}
       {
         if ($2 !~ /^rs[0-9]+$/)
           $2 = $1 ":" $4 ":" $5 ":" $6
         print
       }' "${PREFIX}.bim" > "${PREFIX}.bim.tmp"
  mv "${PREFIX}.bim.tmp" "${PREFIX}.bim"

  awk '{count[$2]++}
       END {
         for (id in count)
           if (count[id] > 1) {
             print "Duplicate variant ID remains: " id > "/dev/stderr"
             exit 1
           }
       }' "${PREFIX}.bim"

  "$PLINK2" --bfile "$PREFIX" \
    --make-grm-bin \
    --pca approx "$NPC" \
    --out "$PREFIX" \
    --threads "$THREADS" \
    --memory "$MEM_MB" >> "$LOG" 2>&1

  awk -v END_COL="$((NPC + 2))" 'BEGIN {OFS=" "}
       NR==1 && $1 ~ /^#?FID$/ {next}
       {
         for (i=1; i<=END_COL; i++)
           printf "%s%s", $i, (i==END_COL ? ORS : OFS)
       }' "${PREFIX}.eigenvec" > "${PREFIX}_adjhe.eigenvec"

  test -s "${PREFIX}_adjhe.eigenvec"
  awk -v END_COL="$((NPC + 2))" 'NF != END_COL {exit 1}' "${PREFIX}_adjhe.eigenvec"

  echo "Cached GRM/PCA written to ${PREFIX}" | tee -a "$LOG"
done
