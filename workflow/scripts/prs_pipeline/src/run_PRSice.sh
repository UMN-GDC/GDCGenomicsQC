#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=30GB
#SBATCH --time=1:00:00
#SBATCH -p msismall
#SBATCH -o PRSice2.out
#SBATCH --job-name PRSice2



path_PRSice=/projects/standard/gdc/public/prs_methods/ref/PRSice/PRSice.R
sum_stats_file=/projects/standard/gdc/public/prs_methods/data/simulated_1000G/gwas/AFR_CT_PRSice2_summary_stat_file.txt
study_sample=/projects/standard/gdc/public/prs_methods/data/simulated_1000G/anc1_plink_files/AFR_simulation_study_sample
phenotype_info_file=/projects/standard/gdc/public/prs_methods/data/simulated_1000G/prs_pipeline/viprs/study_sample_pheno.txt
output_path=/projects/standard/gdc/public/prs_methods/data/simulated_1000G
path_PRSice_linux=/projects/standard/gdc/public/prs_methods/ref/PRSice/PRSice_linux

final_output_dir=${output_path}/prs_pipeline/PRSice

pushd /projects/standard/gdc/public/prs_methods/ref/PRSice
  Rscript ${path_PRSice} \
    --prsice ${path_PRSice_linux} \
    --base ${sum_stats_file} \
    --target ${study_sample} \
    --binary-target F \
    --pheno ${phenotype_info_file} \
    --stat beta \
    --beta \
    --out ${final_output_dir}/PRSice

popd
