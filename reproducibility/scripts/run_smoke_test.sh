#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
run_dir=${1:-"$repo_root/.reproducibility-run"}
mkdir -p "$run_dir"

python3 "$repo_root/reproducibility/scripts/generate_fixture.py" \
  --output-dir "$repo_root/reproducibility/fixtures"

config="$run_dir/smoke-test.yaml"
sed -e "s|@REPO_ROOT@|$repo_root|g" -e "s|@RUN_DIR@|$run_dir|g" \
  "$repo_root/reproducibility/config/smoke-test.yaml.in" > "$config"

cd "$repo_root/workflow"
snakemake --snakefile Snakefile --configfile "$config" --use-conda --cores 2 run_initialQC
python3 "$repo_root/reproducibility/scripts/validate_outputs.py" \
  --output-dir "$run_dir/output" \
  --expectations "$repo_root/reproducibility/expected_outputs.json"
