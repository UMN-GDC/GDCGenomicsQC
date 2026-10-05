# CT-SLEB joint multi-ancestry PRS rule.
#
# Runs the CT-SLEB workflow (2D clumping/thresholding + empirical Bayes +
# super-learning) via run_ctsleb_wrapper.R against the CTSLEB R sources
# vendored at workflow/scripts/CTSLEB/R (or an override set through
# prsMethods.multi_ctsleb.software_dir).
#
# Reference panels ("LD reference" PLINK prefixes) default to the ANC1/ANC2
# study samples with the study sample used both as the clumping LD reference
# and the target test cohort; point prsMethods.multi_ctsleb.target_ref_plink
# / aux_ref_plink at (e.g.) 1000G reference panels for real analyses.

ctsleb_soft_dir = config.get("prsMethods", {}).get("multi_ctsleb", {}).get(
    "software_dir", str(Path(workflow.basedir) / "scripts" / "CTSLEB" / "R"))
ctsleb_anc2_prefix = PRS_OUT_DIR / "anc2_plink_files" / f"{PRS_ANC2}_simulation_study_sample"

rule runMultiAncestryCTSLEB:
    log:
        OUT_DIR / "logs" / "runMultiAncestryCTSLEB.log",
    container:
        "oras://ghcr.io/mainsqu33ze/gdcgenomicsqc/prs:latest"
    conda:
        "../../envs/prs.yml"
    threads: 4
    resources:
        nodes=1,
        mem_mb=32000,
        runtime=720,
    input:
        resources=rules.preparePRSMethodResources.output.ready,
        ss_target=rules.prepareCTSLEBSumstats.output.target_ss,
        ss_training=rules.prepareCTSLEBSumstats.output.training_ss,
        study_bed=rules.preparePRSInputs.output.study_bed,
        study_anc2_bed=rules.preparePRSInputs.output.study_anc2_bed,
        tuning_pheno=rules.prepareCTSLEBPhenotypes.output.tuning,
        validation_pheno=rules.prepareCTSLEBPhenotypes.output.validation,
    output:
        val_r2=PRS_METHOD_RUN_DIR / "multi_ctsleb" / f"{PRS_ANC1}_ctsleb_sl_validation_r2.txt",
        final_betas=PRS_METHOD_RUN_DIR / "multi_ctsleb" / f"{PRS_ANC1}_ctsleb_final_coefficients.txt",
        best_snps=PRS_METHOD_RUN_DIR / "multi_ctsleb" / f"{PRS_ANC1}_ctsleb_best_snps.tsv",
        done=PRS_METHOD_RUN_DIR / "multi_ctsleb.done",
    params:
        out_dir=PRS_METHOD_RUN_DIR / "multi_ctsleb",
        out_prefix=f"{PRS_ANC1}_ctsleb",
        study_prefix=PRS_STUDY_PLINK,
        study_anc2_prefix=ctsleb_anc2_prefix,
        target_ref_plink=lambda wildcards: (
            config.get("prsMethods", {}).get("multi_ctsleb", {}).get("target_ref_plink", "")
            or PRS_STUDY_PLINK
        ),
        aux_ref_plink=lambda wildcards: (
            config.get("prsMethods", {}).get("multi_ctsleb", {}).get("aux_ref_plink", "")
            or ctsleb_anc2_prefix
        ),
        plink19=PRS_PLINK_BIN,
        plink2=PRS_PLINK2_BIN,
        ctsleb_src=ctsleb_soft_dir,
        wrapper=Path(workflow.basedir) / "scripts" / "run_ctsleb_wrapper.R",
    shell:
        """
        set -euo pipefail

        mkdir -p {params.out_dir}

        Rscript {params.wrapper} \
            --target-sumstats {input.ss_target} \
            --aux-sumstats {input.ss_training} \
            --target-ref-plink {params.target_ref_plink} \
            --aux-ref-plink {params.aux_ref_plink} \
            --target-test-plink {params.study_prefix} \
            --plink19 {params.plink19} \
            --plink2 {params.plink2} \
            --results-dir {params.out_dir} \
            --out-prefix {params.out_prefix} \
            --tuning-pheno {input.tuning_pheno} \
            --validation-pheno {input.validation_pheno} \
            --ctsleb-src {params.ctsleb_src} \
            > {log} 2>&1

        touch {output.done}
        """