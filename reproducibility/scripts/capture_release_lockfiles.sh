#!/usr/bin/env bash
# Run after selecting a tagged pipeline release. The resulting files capture
# source, environments, and container identities for reproducible execution.
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
out=${1:-"$repo_root/reproducibility/release-locks"}
mkdir -p "$out/conda" "$out/containers"
git -C "$repo_root" rev-parse HEAD > "$out/source-commit.txt"
git -C "$repo_root" status --porcelain=v1 > "$out/source-status.txt"
date -u +%FT%TZ > "$out/captured-at-utc.txt"

command -v conda-lock >/dev/null || { echo "Install conda-lock before release capture." >&2; exit 2; }
for env in "$repo_root"/envs/*.yml; do
  name=$(basename "${env%.yml}")
  conda-lock lock --kind explicit --file "$env" --filename-template "$out/conda/${name}-{platform}.lock"
done

command -v apptainer >/dev/null || { echo "Install Apptainer before release capture." >&2; exit 2; }
rg -No '(oras|docker)://[^" ]+' "$repo_root/workflow" | cut -d: -f3- | sort -u > "$out/container-uris.txt"
while IFS= read -r image; do
  safe=$(printf '%s' "$image" | tr '/:' '__')
  apptainer pull "$out/containers/${safe}.sif" "$image"
  sha256sum "$out/containers/${safe}.sif"
done < "$out/container-uris.txt" > "$out/container-sha256.txt"
