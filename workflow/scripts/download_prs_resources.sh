#!/usr/bin/env bash
set -euo pipefail

RESOURCE_DIR=""
PRSCSX_REF_DIR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --resource-dir)
      [[ $# -ge 2 ]] || { echo "Missing value for $1" >&2; exit 2; }
      RESOURCE_DIR="$2"; shift 2 ;;
    --prscsx-ref-dir)
      [[ $# -ge 2 ]] || { echo "Missing value for $1" >&2; exit 2; }
      PRSCSX_REF_DIR="$2"; shift 2 ;;
    -h|--help) echo "Usage: $0 --resource-dir DIR [--prscsx-ref-dir DIR]"; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done

[[ -n "$RESOURCE_DIR" ]] || { echo "Missing --resource-dir" >&2; exit 2; }

mkdir -p "$RESOURCE_DIR"/ld{,/prs_cs,/prs_csx,/ldpred2_lassosum2,/ct_sleb,/prosper,/sdprs}

download_if_needed() {
  local url="$1"
  local dest="$2"
  [[ -n "$url" ]] || return 0
  if [[ -e "$dest" ]]; then
    echo "Already exists: $dest"
    return 0
  fi
  mkdir -p "$(dirname "$dest")"
  echo "Downloading $url -> $dest"
  curl -L --fail --retry 3 "$url" -o "$dest"
}

extract_tar_if_present() {
  local archive="$1"
  local dest="$2"

  [[ -f "$archive" ]] || return 0
  if [[ -n "$(find "$dest" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null)" ]]; then
    echo "Already extracted: $dest"
    return 0
  fi

  mkdir -p "$dest"
  echo "Extracting $archive -> $dest"
  tar -xzf "$archive" -C "$dest"
}

link_if_exists() {
  local src="$1"
  local dest="$2"
  [[ -n "$src" && -e "$src" ]] || return 0
  ln -sfn "$src" "$dest"
  echo "Linked $dest -> $src"
}

link_if_exists "$PRSCSX_REF_DIR" "$RESOURCE_DIR/ld/prs_csx/ref"

download_if_needed "${PRSCS_LD_URL:-}" "$RESOURCE_DIR/ld/prs_cs/prs_cs_ld_reference.tar.gz"
download_if_needed "${PRSCSX_LD_URL:-}" "$RESOURCE_DIR/ld/prs_csx/prs_csx_ld_reference.tar.gz"

extract_tar_if_present "$RESOURCE_DIR/ld/prs_cs/prs_cs_ld_reference.tar.gz" "$RESOURCE_DIR/ld/prs_cs/ref"
extract_tar_if_present "$RESOURCE_DIR/ld/prs_csx/prs_csx_ld_reference.tar.gz" "$RESOURCE_DIR/ld/prs_csx/ref"

# Create marker file
touch "$RESOURCE_DIR/resources.ready"
echo "PRS resource layout ready: $RESOURCE_DIR"
