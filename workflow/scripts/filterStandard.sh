#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 INPUT_PFILE_PREFIX OUTPUT_DIRECTORY THREADS" >&2
  exit 2
fi

INPUT=$1
STAGE=$2
THREADS=$3
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
TEMP="$STAGE/intermediates/standard_filter"

mkdir -p "$TEMP"

MAF_PREFIX="$TEMP/maf_filtered"
HWE_PREFIX="$TEMP/hwe_filtered"
HET_LD_PREFIX="$TEMP/heterozygosity_ld"
HET_PREFIX="$STAGE/heterozygosity"
STANDARD_PREFIX="$STAGE/standardFilter"

echo "[standard QC] Input sample and variant counts"
plink2 --pfile "$INPUT" --write-samples --write-snplist \
  --out "$TEMP/input_counts" --threads "$THREADS"

echo "[standard QC] Minor-allele frequency"
plink2 --pfile "$INPUT" --freq --out "$STAGE/MAF_check" \
  --threads "$THREADS"
plink2 --pfile "$INPUT" --maf 0.01 --make-pgen --out "$MAF_PREFIX" \
  --threads "$THREADS"

echo "[standard QC] Hardy-Weinberg equilibrium"
plink2 --pfile "$MAF_PREFIX" --hardy --out "$TEMP/hardy_report" \
  --threads "$THREADS"

# Retain a diagnostic report of variants with HWE P < 1e-5. PLINK2 writes
# .hardy rather than the PLINK1 .hwe filename.
awk 'BEGIN {OFS="\t"}
  NR==1 {
    for (i=1; i<=NF; i++) if ($i=="P") p=i
    if (!p) {
      print "ERROR: P column absent from PLINK2 .hardy report" > "/dev/stderr"
      exit 1
    }
    print
    next
  }
  $p < 1e-5 {print}' "$TEMP/hardy_report.hardy" > "$STAGE/zoomhwe.hardy"

plink2 --pfile "$MAF_PREFIX" --hwe 1e-6 --make-pgen \
  --out "$HWE_PREFIX" --threads "$THREADS"

echo "[standard QC] LD pruning for heterozygosity assessment"
if [[ -s "$SCRIPT_DIR/inversion.txt" ]]; then
  plink2 --pfile "$HWE_PREFIX" \
    --exclude range "$SCRIPT_DIR/inversion.txt" \
    --indep-pairwise 50 5 0.2 \
    --out "$HET_LD_PREFIX" --threads "$THREADS"
else
  plink2 --pfile "$HWE_PREFIX" \
    --indep-pairwise 50 5 0.2 \
    --out "$HET_LD_PREFIX" --threads "$THREADS"
fi

test -s "$HET_LD_PREFIX.prune.in"

plink2 --pfile "$HWE_PREFIX" --extract "$HET_LD_PREFIX.prune.in" \
  --het --out "$HET_PREFIX" --threads "$THREADS"

test -s "$HET_PREFIX.het"

echo "[standard QC] Heterozygosity outliers (mean +/- 3 SD)"
awk 'NR==1 {
    for (i=1; i<=NF; i++) if ($i=="F") f=i
    if (!f) {
      print "ERROR: F column absent from .het report" > "/dev/stderr"
      exit 1
    }
    next
  }
  $f!="NA" {n++; sum+=$f; sumsq+=$f*$f}
  END {
    if (n < 2) exit 1
    mean=sum/n
    variance=(sumsq-(sum*sum/n))/(n-1)
    if (variance < 0 && variance > -1e-15) variance=0
    print mean,sqrt(variance),n
  }' "$HET_PREFIX.het" > "$TEMP/heterozygosity.stats"

read -r HET_MEAN HET_SD HET_N < "$TEMP/heterozygosity.stats"

awk -v mean="$HET_MEAN" -v sd="$HET_SD" 'BEGIN {OFS="\t"}
  NR==1 {
    for (i=1; i<=NF; i++) if ($i=="F") f=i
    print "#FID","IID"
    next
  }
  $f!="NA" && ($f < mean-3*sd || $f > mean+3*sd) {print $1,$2}' \
  "$HET_PREFIX.het" > "$STAGE/het_fail_ind.txt"

N_HET_FAIL=$(( $(wc -l < "$STAGE/het_fail_ind.txt") - 1 ))
echo "Heterozygosity observations: $HET_N"
echo "Heterozygosity mean: $HET_MEAN"
echo "Heterozygosity SD: $HET_SD"
echo "Heterozygosity outliers removed: $N_HET_FAIL"

if [[ "$N_HET_FAIL" -gt 0 ]]; then
  plink2 --pfile "$HWE_PREFIX" --remove "$STAGE/het_fail_ind.txt" \
    --make-pgen --out "$STANDARD_PREFIX" --threads "$THREADS"
else
  plink2 --pfile "$HWE_PREFIX" --make-pgen \
    --out "$STANDARD_PREFIX" --threads "$THREADS"
fi

echo "[standard QC] Final marker list and LD-pruned dataset"
plink2 --pfile "$STANDARD_PREFIX" --write-snplist \
  --out "$STANDARD_PREFIX" --threads "$THREADS"
plink2 --pfile "$STANDARD_PREFIX" --indep-pairwise 500 10 0.1 \
  --out "$STANDARD_PREFIX" --threads "$THREADS"

test -s "$STANDARD_PREFIX.prune.in"

plink2 --pfile "$STANDARD_PREFIX" --extract "$STANDARD_PREFIX.prune.in" \
  --make-pgen --out "$STANDARD_PREFIX.LDpruned" --threads "$THREADS"

echo "[standard QC] Complete"
