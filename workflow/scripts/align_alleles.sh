#!/bin/bash
# =============================================================================
# align_alleles.sh — Compare study variant alleles to reference panel by
# position, identify strand flips, and write a flip list for PLINK2 --flip.
#
# Usage:
#   align_alleles.sh <study.pvar> <ref.pvar> <flip_list.txt> <report.txt>
#
# Input .pvar files are tab-separated with columns:
#   #CHROM  POS     ID      REF     ALT
#
# Matching is by (chromosome, position). Allele sets are compared:
#   - exact match: same REF/ALT (any order) → pass
#   - strand flip: complemented alleles → write study ID to flip list
#   - mismatch: different alleles → logged as excluded
# =============================================================================

set -euo pipefail

STUDY_PVAR=$1
REF_PVAR=$2
FLIP_OUT=$3
REPORT_OUT=$4

# Safety: clear output files before appending
> "$FLIP_OUT"
> "$REPORT_OUT"

# Count total ref variants (exclude header line)
TOTAL_REF=$(awk '!/^#/{count++} END{print count+0}' "$REF_PVAR")

awk -v flip_out="$FLIP_OUT" -v report_out="$REPORT_OUT" -v total_ref="$TOTAL_REF" '
BEGIN {
    OFS = "\t"
    comp["A"] = "T"; comp["T"] = "A"; comp["C"] = "G"; comp["G"] = "C"
    comp["a"] = "t"; comp["t"] = "a"; comp["c"] = "g"; comp["g"] = "c"
}
NR==FNR {
    if (/^#/) next
    ref[$1,$2] = $4 OFS $5
    next
}
!/^#/ {
    key = $1 OFS $2
    if (key in ref) {
        split(ref[key], r, OFS)
        s_ref = toupper($4); s_alt = toupper($5)
        r_ref = toupper(r[1]); r_alt = toupper(r[2])

        # Same allele set (any REF/ALT ordering)
        if ((s_ref == r_ref && s_alt == r_alt) || (s_ref == r_alt && s_alt == r_ref)) {
            exact++
        }
        # Complement match — strand flip needed
        else if ((s_ref == comp[r_ref] && s_alt == comp[r_alt]) || (s_ref == comp[r_alt] && s_alt == comp[r_ref])) {
            flip++
            print $3 > flip_out   # variant ID for --flip
        }
        # Genuine allele mismatch
        else {
            mismatch++
        }
    } else {
        study_only++
    }
}
END {
    print "exact_match:", exact+0 > report_out
    print "strand_flip:", flip+0 > report_out
    print "allele_mismatch:", mismatch+0 > report_out
    print "study_only:", study_only+0 > report_out
    printf "ref_only: %d\n", total_ref - exact - flip - mismatch > report_out
}
' "$REF_PVAR" "$STUDY_PVAR"

N_FLIP=$(wc -l < "$FLIP_OUT" 2>/dev/null || echo 0)

# Read stats from report for stdout
EXACT=$(grep "^exact_match:" "$REPORT_OUT" | cut -d' ' -f2)
FLIP=$(grep "^strand_flip:" "$REPORT_OUT" | cut -d' ' -f2)
MISMATCH=$(grep "^allele_mismatch:" "$REPORT_OUT" | cut -d' ' -f2)
STUDY_ONLY=$(grep "^study_only:" "$REPORT_OUT" | cut -d' ' -f2)
REF_ONLY=$(grep "^ref_only:" "$REPORT_OUT" | cut -d' ' -f2)

echo "[align_alleles] matched $EXACT exact, $FLIP flipped, $MISMATCH mismatched, $STUDY_ONLY study-only, $REF_ONLY ref-only"
echo "[align_alleles] wrote $N_FLIP variants to flip list"
