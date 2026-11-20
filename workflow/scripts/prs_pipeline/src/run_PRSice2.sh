#!/bin/bash
summary_stats_file=$1
study_sample=$2
binary_flag=$3
phenotype_info_file=$4
output_path=$5
path_repo=$6
output_prefix=$7   # optional: overrides the default output prefix

if [[ -n "${output_prefix:-}" ]]; then
    final_output="${output_prefix}"
else
    final_output="${output_path}/prs_pipeline/PRSice2/prs_method"
fi
mkdir -p "$final_output"

# run PRSice-2 on simulated data--note that this is the temporary data and is unlikely to produce correct results (use this more for formatting purposes)
PRSICE_CMD=${PRSICE_CMD:-PRSice}
PLINK_CMD=${PLINK_CMD:-plink}
${PRSICE_CMD} \
    --base ${summary_stats_file} \
    --target ${study_sample} \
    --binary-target ${binary_flag} \
    --pheno ${phenotype_info_file} \
    --stat beta \
    --beta \
    --out ${final_output}/PRSice2_outputs

# Generate clumped SNP list for test evaluation (used by score_test.sh)
# PLINK clumping with same default parameters as PRSice2 / C+T (r2=0.1, kb=250)
${PLINK_CMD} \
    --bfile "${study_sample}" \
    --clump-p1 1 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --clump "${summary_stats_file}" \
    --clump-snp-field SNP \
    --clump-field P \
    --allow-no-sex \
    --out "${final_output}/PRSice2_clump" 2>/dev/null

awk 'NR!=1{print $3}' "${final_output}/PRSice2_clump.clumped" > "${final_output}/PRSice2_outputs.snps" 2>/dev/null || true
