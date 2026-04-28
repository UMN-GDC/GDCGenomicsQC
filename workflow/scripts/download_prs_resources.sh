#!/usr/bin/env bash
set -euo pipefail

RESOURCE_DIR=""
PRSCSX_REF_DIR=""
PLINK2=""
DOWNLOAD_SOFTWARE="false"

usage() {
  cat <<'USAGE'
Usage:
  download_prs_resources.sh --resource-dir DIR [--prscsx-ref-dir DIR] [--plink2 PATH] [--download-software]

Creates a standard PRS resource layout. The script prefers symlinks to existing
MSI/shared resources and leaves method-specific downloads configurable because
several PRS reference panels have separate licenses or large external archives.

Optional environment variables:
  PRSCS_LD_URL       URL to PRS-CS LD reference archive.
  PRSCSX_LD_URL      URL to PRS-CSx LD reference archive.
  PRSICE_URL         URL to PRSice-2 archive.
  CTSLEB_URL         URL to CT-SLEB software archive/repository tarball.
  PROSPER_URL        URL to PROSPER software archive/repository tarball.
  SDPRS_URL          URL to SDPRS/SDPRX software archive/repository tarball.
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --resource-dir) RESOURCE_DIR="$2"; shift 2 ;;
    --prscsx-ref-dir) PRSCSX_REF_DIR="$2"; shift 2 ;;
    --plink2) PLINK2="$2"; shift 2 ;;
    --download-software) DOWNLOAD_SOFTWARE="true"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -n "$RESOURCE_DIR" ]] || { usage >&2; exit 2; }

mkdir -p "$RESOURCE_DIR"/{software,ld,logs}
mkdir -p "$RESOURCE_DIR"/ld/{prs_cs,prs_csx,ldpred2_lassosum2,ct_sleb,prosper,sdprs}

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

clone_if_needed() {
  local repo="$1"
  local dest="$2"

  if [[ -d "$dest/.git" ]]; then
    echo "Already cloned: $dest"
    return 0
  fi

  command -v git >/dev/null 2>&1 || {
    echo "git is required for --download-software" >&2
    exit 127
  }

  mkdir -p "$(dirname "$dest")"
  echo "Cloning $repo -> $dest"
  git clone --depth 1 "$repo" "$dest"
}

link_if_exists() {
  local src="$1"
  local dest="$2"
  [[ -n "$src" && -e "$src" ]] || return 0
  ln -sfn "$src" "$dest"
  echo "Linked $dest -> $src"
}

link_if_exists "$PRSCSX_REF_DIR" "$RESOURCE_DIR/ld/prs_csx/ref"
link_if_exists "$PLINK2" "$RESOURCE_DIR/software/plink2"

if [[ "$DOWNLOAD_SOFTWARE" == "true" ]]; then
  clone_if_needed "https://github.com/getian107/PRScs.git" "$RESOURCE_DIR/software/PRScs"
  clone_if_needed "https://github.com/getian107/PRScsx.git" "$RESOURCE_DIR/software/PRScsx"
  clone_if_needed "https://github.com/andrewhaoyu/CTSLEB.git" "$RESOURCE_DIR/software/CTSLEB"
  clone_if_needed "https://github.com/Jingning-Zhang/PROSPER.git" "$RESOURCE_DIR/software/PROSPER"
  clone_if_needed "https://github.com/eldronzhou/SDPRX.git" "$RESOURCE_DIR/software/SDPRX"
fi

download_if_needed "${PRSCS_LD_URL:-}" "$RESOURCE_DIR/ld/prs_cs/prs_cs_ld_reference.tar.gz"
download_if_needed "${PRSCSX_LD_URL:-}" "$RESOURCE_DIR/ld/prs_csx/prs_csx_ld_reference.tar.gz"
download_if_needed "${PRSICE_URL:-}" "$RESOURCE_DIR/software/prsice.tar.gz"
download_if_needed "${CTSLEB_URL:-}" "$RESOURCE_DIR/software/ctsleb.tar.gz"
download_if_needed "${PROSPER_URL:-}" "$RESOURCE_DIR/software/prosper.tar.gz"
download_if_needed "${SDPRS_URL:-}" "$RESOURCE_DIR/software/sdprs.tar.gz"

extract_tar_if_present "$RESOURCE_DIR/ld/prs_cs/prs_cs_ld_reference.tar.gz" "$RESOURCE_DIR/ld/prs_cs/ref"
extract_tar_if_present "$RESOURCE_DIR/ld/prs_csx/prs_csx_ld_reference.tar.gz" "$RESOURCE_DIR/ld/prs_csx/ref"

cat > "$RESOURCE_DIR/resources.ready" <<EOF
resource_dir="$RESOURCE_DIR"
created_at="$(date -Iseconds)"
prscsx_ref_dir="$PRSCSX_REF_DIR"
plink2="$PLINK2"
download_software="$DOWNLOAD_SOFTWARE"
EOF

echo "PRS resource layout ready: $RESOURCE_DIR"
