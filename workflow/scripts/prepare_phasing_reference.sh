#!/usr/bin/env bash
# Rename the contigs of a VCF to the canonical bare form ("20" instead of
# "chr20").
#
# The 1KG high-coverage panel ships "chr"-prefixed contigs, while the study
# PLINK/PGEN data and the genetic maps use bare contigs. Every consumer in the
# local-ancestry chain (site intersection, shapeit4 --region, RFMIX/gnomix
# reference) must agree on one namespace, so the rename happens here, once,
# from the contig names in the file's own header.
#
# Usage: prepare_phasing_reference.sh <in.vcf.gz> <out.vcf.gz> [--index]
set -euo pipefail

IN=$1
OUT=$2
INDEX=${3:-}

if [[ -z "${IN:-}" || -z "${OUT:-}" ]]; then
    echo "Usage: $0 <in.vcf.gz> <out.vcf.gz> [--index]" >&2
    exit 2
fi

MAP="${OUT}.rename.txt"

# Contig names come from the header only (no full-file scan).
bcftools view -h "$IN" \
    | sed -n 's/^##contig=<ID=\([^,]*\).*/\1/p' \
    | sort -u \
    | awk 'BEGIN {OFS="\t"} {c=$0; sub(/^chr/, "", c); if (c != $0) print $0, c}' > "$MAP"

if [[ -s "$MAP" ]]; then
    bcftools annotate --rename-chrs "$MAP" "$IN" -Oz -o "$OUT"
else
    # Already bare: keep the file as it is (avoids a pointless re-encode).
    bcftools view -Oz -o "$OUT" "$IN"
fi
rm -f "$MAP"

if [[ "$INDEX" == "--index" ]]; then
    bcftools index -f "$OUT"
    bcftools index -t -f "$OUT"
fi
