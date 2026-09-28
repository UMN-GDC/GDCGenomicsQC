#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=30GB
#SBATCH --time=1:00:00
#SBATCH -p msismall
#SBATCH -o CT_prs.out
#SBATCH --job-name CT_prs

# prior to running, make sure that correct summary statistics file has been generated; look into whether or not the file needs to be stored in the same place as the data or if you can read it in separately
## Inputs
# plink bfiles from study population of interest
# GWAS summary stats file 
# Output path

load_config() {
    local file_path="$1"
    if [[ -f "$file_path" ]]; then
        echo "Loading configuration from: $file_path"
        # Source the file. Variables set in the config will override defaults.
        # Note: 'source' is used for shell-readable KEY="VALUE" files.
        source "$file_path"
        return 0
    else
        echo "ERROR: Configuration file not found at $file_path" >&2
        return 1
    fi
}

## Defaults
study_sample=EUR_simulation_study_sample
sum_stats_file=sumstats_CT_PRSice2.txt
phenotype_info_file=phenotype_info.txt
gwas_pca_eigenvec_file=""
output_path=/projects/standard/gdc/public/data/simulated_1000G
path_prs_pipeline=/projects/standard/gdc/public/prs_methods/prs_pipeline

#### # ---- Parse args ----
# Check if the number of arguments is 0
if [ "$#" -eq 0 ]; then
    echo "Error: No config file provided."
    echo "Usage: $0 <config>"
    exit 1
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    --c) config_file="$2"; shift 2; break ;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown option: $1"; usage; exit 1;;
  esac
done

if [[ -n "$config_file" ]]; then
    # Case 1: --config was specified. Load it, and it overrides ALL defaults.
    load_config "$config_file" || exit 1
else
    # Case 2: No --config was specified.
    # The variables set by the individual flags in step 3 are used.
    echo "No --c file specified. Using command-line arguments and defaults."
fi


final_output_dir=${final_output_dir:-${output_path}/CT}

# source /projects/standard/gdc/public/envs/load_miniconda3.sh

######### Script start ###########
mkdir -p ${final_output_dir}/temp

pushd ${final_output_dir}/temp

# plink commands
plink \
    --bfile ${study_sample} \
    --clump-p1 1 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --clump ${sum_stats_file} \
    --clump-snp-field SNP \
    --clump-field P \
    --allow-no-sex \
    --out temp

awk 'NR!=1{print $3}' temp.clumped >  temp.valid.snp
awk '{print $3,$8}' ${sum_stats_file}  > SNP.pvalue
echo "0.001 0 0.001" > range_list
echo "0.05 0 0.05" >> range_list
echo "0.1 0 0.1" >> range_list
echo "0.2 0 0.2" >> range_list
echo "0.3 0 0.3" >> range_list
echo "0.4 0 0.4" >> range_list
echo "0.5 0 0.5" >> range_list

awk 'NR==1 {print "SNP pvalue"; next} {print $1, $8}' OFS="\t" ${sum_stats_file} > SNP.pvalue

plink \
    --bfile ${study_sample} \
    --score ${sum_stats_file} 1 4 6 header \
    --q-score-range range_list SNP.pvalue \
    --extract temp.valid.snp \
    --allow-no-sex \
    --out temp
plink \
    --bfile ${study_sample} \
    --indep-pairwise 200 50 0.25 \
    --allow-no-sex \
    --out temp

plink \
    --bfile ${study_sample} \
    --extract temp.prune.in \
    --pca 6 \
    --allow-no-sex \
    --out temp

# Use study-sample PCA by default; user-provided gwas_pca_eigenvec_file overrides
pca_file="${gwas_pca_eigenvec_file:-temp.eigenvec}"
if [[ -z "$gwas_pca_eigenvec_file" ]]; then
    pca_file="temp.eigenvec"
    echo "Using study-sample PCA (temp.eigenvec)"
else
    pca_file="$gwas_pca_eigenvec_file"
    echo "Using user-provided PCA file: $pca_file"
fi

# R commands--helps find the p-value threshold that leads to the PRS with the best fit under the clumping and thresholding model
## note that this is approximate (it cannot usually be determined with certainty) but can be approximated by performing a regression 
## between the calculated PRS at a given range of p-values; can then select the PRS that explains the highest phenotypic variance
Rscript ${path_prs_pipeline}/src/CT_script.R \
    ${phenotype_info_file} \
    ${pca_file} \
    temp \
    ${final_output_dir}

