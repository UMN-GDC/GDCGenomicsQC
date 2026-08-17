rule runMultiAncestryCTSLEB:
    input:
        resources=rules.preparePRSMethodResources.output,
        prs_env=OUT_DIR / "prs_inputs_hm3" / "AFR_EUR" / "prs_inputs.env",
        target_ss=OUT_DIR / "prs_inputs_hm3" / "AFR_EUR" / "gwas" / "target_sumstats.txt",
        training_ss=OUT_DIR / "prs_inputs_hm3" / "AFR_EUR" / "gwas" / "training_sumstats.txt",
        target_study=OUT_DIR / "prs_inputs_hm3" / "AFR_EUR" / "anc1_plink_files" / "AFR_simulation_study_sample.bed",
        training_study=OUT_DIR / "prs_inputs_hm3" / "AFR_EUR" / "anc2_plink_files" / "EUR_simulation_study_sample.bed",
    output:
        done=OUT_DIR / "prs_inputs_hm3" / "AFR_EUR" / "method_runs" / "multi_ctsleb.done",
    log:
        OUT_DIR / "logs" / "runMultiAncestryCTSLEB.log",
    params:
        out_dir=OUT_DIR / "prs_inputs_hm3" / "AFR_EUR" / "method_runs" / "multi_ctsleb",
        resource_dir=config.get("prsMethods", {}).get("resource_dir", ""),
        software_dir=lambda wildcards: config.get("prsMethods", {}).get("multi_ctsleb", {}).get("software_dir", ""),
        ld_ref_dir=lambda wildcards: config.get("prsMethods", {}).get("multi_ctsleb", {}).get("ld_ref_dir", ""),
        wrapper=Path(workflow.basedir) / "scripts" / "run_ctsleb_wrapper.R",
        plink19=config.get("prsMethods", {}).get("multi_ctsleb", {}).get("plink19", "/projects/standard/gdc/public/plink"),
        plink2=config.get("prsMethods", {}).get("plink2", "/projects/standard/gdc/public/plink2"),
        target_ref_plink=lambda wildcards: config.get("prsMethods", {}).get("multi_ctsleb", {}).get("target_ref_plink", ""),
        aux_ref_plink=lambda wildcards: config.get("prsMethods", {}).get("multi_ctsleb", {}).get("aux_ref_plink", ""),
        tuning_pheno=lambda wildcards: config.get("prsMethods", {}).get("multi_ctsleb", {}).get("tuning_pheno", ""),
        validation_pheno=lambda wildcards: config.get("prsMethods", {}).get("multi_ctsleb", {}).get("validation_pheno", ""),
        out_prefix=lambda wildcards: config.get("prsMethods", {}).get("multi_ctsleb", {}).get("out_prefix", "AFR_EUR_ctsleb"),
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.out_dir}

        Rscript {params.wrapper} \
          --target-sumstats {input.target_ss} \
          --aux-sumstats {input.training_ss} \
          --target-ref-plink {params.target_ref_plink} \
          --aux-ref-plink {params.aux_ref_plink} \
          --target-test-plink {input.target_study[:-4]} \
          --plink19 {params.plink19} \
          --plink2 {params.plink2} \
          --results-dir {params.out_dir} \
          --out-prefix {params.out_prefix} \
          --tuning-pheno {params.tuning_pheno} \
          --validation-pheno {params.validation_pheno} \
          > {log} 2>&1

        touch {output.done}
        """
