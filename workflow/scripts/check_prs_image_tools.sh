#!/usr/bin/env bash
# Verify that the tools the PRS method rules need are provided by the
# oras://ghcr.io/mainsqu33ze/gdcgenomicsqc/prs:latest image (PATH + baked
# /opt/gdcgenomicsqc software). Replaces the old download_prs_resources.sh:
# software is no longer cloned onto the host. PROSPER is the documented
# exception (singleprshelper image + host software_dir/plink2 overrides) and
# is intentionally not checked here.
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: check_prs_image_tools.sh --out FILE [--prscsx-code DIR] [--ctsleb-r DIR]
USAGE
}

OUT=""
PRSCSX_CODE="/opt/gdcgenomicsqc/software/PRScsx"
CTSLEB_R="/opt/gdcgenomicsqc/CTSLEB/R"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --out) OUT="$2"; shift 2 ;;
    --prscsx-code) PRSCSX_CODE="$2"; shift 2 ;;
    --ctsleb-r) CTSLEB_R="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ -z "$OUT" ]]; then
  usage >&2
  exit 2
fi

fail() { echo "ERROR: $*" >&2; exit 1; }

for tool in PRSice plink plink2 python Rscript; do
  resolved=$(command -v "$tool" || true)
  [[ -n "$resolved" ]] || fail "$tool not found on PATH (expected from the prs image)"
  echo "$tool=$resolved"
done

[[ -f "$PRSCSX_CODE/PRScsx.py" ]] || \
  fail "PRScsx.py not found at $PRSCSX_CODE/PRScsx.py"
echo "prscsx_code=$PRSCSX_CODE"

[[ -f "$CTSLEB_R/helper_CombinePRS.R" ]] || \
  fail "CTSLEB R sources not found at $CTSLEB_R (expected helper_CombinePRS.R)"
echo "ctsleb_r=$CTSLEB_R"

mkdir -p "$(dirname "$OUT")"
{
  echo "image=oras://ghcr.io/mainsqu33ze/gdcgenomicsqc/prs:latest"
  echo "PRSice=$(command -v PRSice)"
  echo "plink=$(command -v plink)"
  echo "plink2=$(command -v plink2)"
  echo "python=$(command -v python)"
  echo "Rscript=$(command -v Rscript)"
  echo "prscsx_code=$PRSCSX_CODE"
  echo "ctsleb_r=$CTSLEB_R"
  echo "prosper=external (singleprshelper image + host software_dir/plink2 overrides)"
  echo "checked_at=$(date -Iseconds)"
} > "$OUT"

echo "OK: PRS image tools verified; marker written to $OUT"
