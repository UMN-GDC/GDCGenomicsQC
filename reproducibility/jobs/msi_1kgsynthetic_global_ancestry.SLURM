#!/usr/bin/env bash
#SBATCH --job-name=gdc-1kg-ancestry
#SBATCH --output=/scratch.global/%u/gdc-1kgsynthetic-validation/master-%j.log
#SBATCH --account=gdc
#SBATCH --time=72:00:00
#SBATCH --mem=4G

set -euo pipefail

# Slurm copies a batch script to its spool directory, so BASH_SOURCE does not
# point to the checked-out repository at execution time. Submit this script
# from the repository root; SLURM_SUBMIT_DIR preserves that directory.
repo_root="${SLURM_SUBMIT_DIR:?Submit this job from the GDCGenomicsQC repository root.}"
if [[ ! -f "$repo_root/workflow/Snakefile" ]]; then
  echo "SLURM_SUBMIT_DIR is not a GDCGenomicsQC checkout: $repo_root" >&2
  exit 2
fi
run_dir="/scratch.global/${USER}/gdc-1kgsynthetic-validation"
mkdir -p "$run_dir"

# MSI interactive logins may initialize an older Conda function. Use the
# current Miniforge executable for Snakemake's environment management.
module load miniforge/26.7
unset -f conda __conda_activate __conda_exe __conda_reactivate __conda_hashr 2>/dev/null || true
hash -r
eval "$(conda shell.bash hook)"
conda activate /projects/standard/gdc/public/envs/snakemake

config="$run_dir/1kgsynthetic-global-ancestry.yaml"
sed "s|@RUN_DIR@|$run_dir|g" \
  "$repo_root/reproducibility/config/msi-1kgsynthetic-global-ancestry.yaml.in" > "$config"

git -C "$repo_root" rev-parse HEAD >"$run_dir/git-commit.txt"
conda --version >"$run_dir/conda-version.txt"
snakemake --version >"$run_dir/snakemake-version.txt"
sha256sum "$config" >"$run_dir/config.sha256"

cd "$repo_root/workflow"
snakemake --profile ../profiles/hpc --configfile "$config" --rerun-incomplete run_classifyAncestry

python3 "$repo_root/reproducibility/scripts/evaluate_global_ancestry.py" \
  --truth /projects/standard/gdc/public/Ref/toyData/1kgSynthetic.popu \
  --classifications "$run_dir/output/01-globalAncestry/ancestry_classifications.tsv" \
  --output-dir "$run_dir/evaluation" \
  --model pca

find "$run_dir/output" -type f -print0 | sort -z | xargs -0 sha256sum >"$run_dir/output.sha256"
