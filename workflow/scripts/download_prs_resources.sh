#!/usr/bin/env bash
set -euo pipefail

RESOURCE_DIR=""
PRSCSX_REF_DIR=""
PLINK2=""
PRS_PIPELINE_DIR=""
PRS_PIPELINE_REF=""
PRS_PIPELINE_SIF=""
PRS_HELPER_SIF=""
DOWNLOAD_SOFTWARE="false"

usage() {
  cat <<'USAGE'
Usage:
  download_prs_resources.sh --resource-dir DIR [--prscsx-ref-dir DIR] [--plink2 PATH] [--prs-pipeline-dir DIR] [--prs-pipeline-ref REF] [--prs-pipeline-sif FILE] [--prs-helper-sif FILE] [--download-software]

Creates a standard PRS resource layout. The script prefers symlinks to existing
MSI/shared resources and leaves method-specific downloads configurable because
several PRS reference panels have separate licenses or large external archives.

--prs-pipeline-dir/--prs-pipeline-ref: resolve and record the absorbed prs_pipeline
method engine (default: the vendored workflow/scripts/prs_pipeline copy; see
VENDORED.md). Must contain run_single_ancestry_PRS_pipeline.sh and src/. The ref
(branch sandbox_multi_pheno @ f71cf4f) is recorded as provenance; HEAD is
verified against it only when --prs-pipeline-dir is a genuine upstream
prs_pipeline git checkout (origin remote contains prs_pipeline), since a vendored
copy's git HEAD is GDCGenomicsQC's, not upstream's.

--prs-pipeline-sif/--prs-helper-sif: SIF container paths (recorded for phase 7,
not yet consumed). The PRS compute rules run inside the GDCGenomicsQC-owned
`prs` image (built from envs/prs.def, bakes the engine + pinned method repos).
Pull with:
  apptainer pull oras://ghcr.io/mainsqu33ze/gdcgenomicsqc/prs:latest
  apptainer pull oras://ghcr.io/mainsqu33ze/gdcgenomicsqc/singleprshelper:latest

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
    --resource-dir)
      if [[ $# -ge 2 && -n "$2" && "$2" != --* ]]; then RESOURCE_DIR="$2"; shift 2; else shift 1; fi ;;
    --prscsx-ref-dir)
      if [[ $# -ge 2 && -n "$2" && "$2" != --* ]]; then PRSCSX_REF_DIR="$2"; shift 2; else shift 1; fi ;;
    --plink2)
      if [[ $# -ge 2 && -n "$2" && "$2" != --* ]]; then PLINK2="$2"; shift 2; else shift 1; fi ;;
    --prs-pipeline-dir)
      if [[ $# -ge 2 && -n "$2" && "$2" != --* ]]; then PRS_PIPELINE_DIR="$2"; shift 2; else shift 1; fi ;;
    --prs-pipeline-ref)
      if [[ $# -ge 2 && -n "$2" && "$2" != --* ]]; then PRS_PIPELINE_REF="$2"; shift 2; else shift 1; fi ;;
    --prs-pipeline-sif)
      if [[ $# -ge 2 && -n "$2" && "$2" != --* ]]; then PRS_PIPELINE_SIF="$2"; shift 2; else shift 1; fi ;;
    --prs-helper-sif)
      if [[ $# -ge 2 && -n "$2" && "$2" != --* ]]; then PRS_HELPER_SIF="$2"; shift 2; else shift 1; fi ;;
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

verify_prs_pipeline() {
  local repo="$1"
  local ref="$2"
  local origin

  [[ -n "$repo" ]] || return 0

  if [[ ! -d "$repo" ]]; then
    echo "ERROR: prs_pipeline tree not found at $repo (check prs_pipeline_path)" >&2
    exit 1
  fi
  for f in run_single_ancestry_PRS_pipeline.sh src/prepare_sumstats.R src/run_CT.sh; do
    if [[ ! -f "$repo/$f" ]]; then
      echo "ERROR: missing vendored $repo/$f" >&2
      exit 1
    fi
  done

  # A vendored copy lives inside the GDCGenomicsQC repo (origin != prs_pipeline),
  # so git HEAD there refers to GDCQC, not the upstream project. Only do upstream
  # pin verification when the configured path is a genuine prs_pipeline checkout.
  origin=$(git -C "$repo" config --get remote.origin.url 2>/dev/null || true)
  if [[ -d "$repo/.git" && "$origin" == *prs_pipeline* ]]; then
    echo "prs_pipeline HEAD: $(git -C "$repo" rev-parse HEAD) ($(git -C "$repo" rev-parse --abbrev-ref HEAD))"
    if [[ -n "$ref" ]]; then
      if git -C "$repo" rev-parse --verify "$ref"^{commit} >/dev/null 2>&1; then
        pinned=$(git -C "$repo" rev-parse "$ref"^{commit})
        head=$(git -C "$repo" rev-parse HEAD)
        if [[ "$head" == "$pinned" ]]; then
          echo "prs_pipeline pinned ref OK ($ref)"
        else
          echo "WARNING: prs_pipeline HEAD ($head) != pinned ref ($ref)" >&2
          echo "  Fix with: git -C $repo checkout $ref" >&2
        fi
      else
        echo "WARNING: pinned ref '$ref' not present locally. Run: git -C $repo fetch --all" >&2
      fi
    fi
  else
    echo "prs_pipeline: vendored copy in GDCGenomicsQC (ref $ref recorded as provenance)"
  fi
}

verify_prs_pipeline "$PRS_PIPELINE_DIR" "$PRS_PIPELINE_REF"

if git -C "$PRS_PIPELINE_DIR" config --get remote.origin.url 2>/dev/null | grep -q prs_pipeline; then
  prs_pipeline_head="$(git -C "$PRS_PIPELINE_DIR" rev-parse HEAD 2>/dev/null || echo unknown)"
  prs_pipeline_vendored="false"
else
  prs_pipeline_head=""
  prs_pipeline_vendored="true"
fi

cat > "$RESOURCE_DIR/resources.ready" <<EOF
resource_dir="$RESOURCE_DIR"
created_at="$(date -Iseconds)"
prscsx_ref_dir="$PRSCSX_REF_DIR"
plink2="$PLINK2"
download_software="$DOWNLOAD_SOFTWARE"
prs_pipeline_dir="$PRS_PIPELINE_DIR"
prs_pipeline_ref="$PRS_PIPELINE_REF"
prs_pipeline_head="$prs_pipeline_head"
prs_pipeline_vendored="$prs_pipeline_vendored"
prs_pipeline_sif="$PRS_PIPELINE_SIF"
prs_helper_sif="$PRS_HELPER_SIF"
EOF

echo "PRS resource layout ready: $RESOURCE_DIR"
