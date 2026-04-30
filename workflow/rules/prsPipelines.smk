PRS_METHODS_CONFIG = config.get("prsMethods", {})
PRS_RESOURCE_DIR = Path(
    PRS_METHODS_CONFIG.get(
        "resource_dir",
        str(OUT_DIR.parent / "prs_resources"),
    )
)
PRS_METHOD_RUN_DIR = PRS_OUT_DIR / "method_runs"

def prs_method_enabled(method):
    """Check if a PRS method is enabled (true/false in config)."""
    return PRS_METHODS_CONFIG.get(method, {}).get("enabled", False)

def get_enabled_prs_methods():
    """Return list of enabled PRS method names."""
    methods = [
        "single_ct", "single_prsice", "single_prscs", "single_ldpred2",
        "single_lassosum2", "multi_ctsleb", "multi_prscsx", "multi_ldpred2",
        "multi_prosper", "multi_sdprs",
    ]
    return [m for m in methods if prs_method_enabled(m)]

def get_method_extra_args(method):
    """Get extra arguments for a PRS method from config."""
    method_config = PRS_METHODS_CONFIG.get(method, {})
    args = []
    for key, flag in (
        ("ld_ref_dir", "--ld-ref-dir"),
        ("ld_ref_prefix", "--ld-ref-prefix"),
        ("ld_matrix_dir", "--ld-matrix-dir"),
    ):
        value = method_config.get(key)
        if value:
            args.extend([flag, str(value)])
    return " ".join(args)


rule preparePRSMethodResources:
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/prs:latest"
    log:
        OUT_DIR / "logs" / "preparePRSMethodResources.log",
    threads: 1
    resources:
        nodes=1,
        mem_mb=4000,
        runtime=60,
    output:
        ready=PRS_RESOURCE_DIR / "resources.ready",
        ld_dir=directory(PRS_RESOURCE_DIR / "ld"),
    params:
        resource_dir=PRS_RESOURCE_DIR,
    shell:
        """
        set -euo pipefail
        mkdir -p "{params.resource_dir}"/ld/prs_cs
        mkdir -p "{params.resource_dir}"/ld/prs_csx
        mkdir -p "{params.resource_dir}"/ld/ldpred2_lassosum2
        mkdir -p "{params.resource_dir}"/ld/ct_sleb
        mkdir -p "{params.resource_dir}"/ld/prosper
        mkdir -p "{params.resource_dir}"/ld/sdprs
        touch "{output.ready}"
        """


rule runSingleAncestryCT:
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/prs:latest"
    log:
        OUT_DIR / "logs" / "runSingleAncestryCT.log",
    threads: 4
    resources:
        nodes=1,
        mem_mb=16000,
        runtime=240,
    input:
        resources=rules.preparePRSMethodResources.output.ready,
        env=rules.preparePRSInputs.output.env,
    output:
        done=PRS_METHOD_RUN_DIR / "single_ct.done",
    params:
        enabled=prs_method_enabled("single_ct"),
        out_dir=PRS_METHOD_RUN_DIR / "single_ct",
        script=Path(workflow.basedir) / "scripts" / "ctsleb.R",
        extra=get_method_extra_args("single_ct"),
    shell:
        """
        if [[ "{params.enabled}" != "True" ]]; then
            echo "Method single_ct not enabled, skipping..."
            touch {output.done}
            exit 0
        fi
        source {input.env}
        mkdir -p {params.out_dir}
        Rscript {params.script} \
            --ss "$target_sumstats_file" \
            --bed "$study_sample_plink.bed" \
            --bim "$study_sample_plink.bim" \
            --fam "$study_sample_plink.fam" \
            --pheno "$target_study_pheno_file" \
            --out-dir {params.out_dir} \
            {params.extra} \
            > {log} 2>&1
        touch {output.done}
        """

rule runSingleAncestryLassosum2:
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/prs:latest"
    log:
        OUT_DIR / "logs" / "runSingleAncestryLassosum2.log",
    threads: 4
    resources:
        nodes=1,
        mem_mb=64000,
        runtime=720,
    input:
        resources=rules.preparePRSMethodResources.output.ready,
        env=rules.preparePRSInputs.output.env,
    output:
        done=PRS_METHOD_RUN_DIR / "single_lassosum2.done",
    params:
        enabled=prs_method_enabled("single_lassosum2"),
        out_dir=PRS_METHOD_RUN_DIR / "single_lassosum2",
        script=Path(workflow.basedir) / "scripts" / "lassosum2.R",
        extra=get_method_extra_args("single_lassosum2"),
    shell:
        """
        if [[ "{params.enabled}" != "True" ]]; then
            echo "Method single_lassosum2 not enabled, skipping..."
            touch {output.done}
            exit 0
        fi
        source {input.env}
        mkdir -p {params.out_dir}
        Rscript {params.script} \
            --sumstats "$target_sumstats_file" \
            --bed "$study_sample_plink.bed" \
            --bim "$study_sample_plink.bim" \
            --fam "$study_sample_plink.fam" \
            --pheno "$target_study_pheno_file" \
            --out-dir {params.out_dir} \
            {params.extra} \
            > {log} 2>&1
        touch {output.done}
        """


rule runSingleAncestryPRSice:
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/prs:latest"
    log:
        OUT_DIR / "logs" / "runSingleAncestryPRSice.log",
    threads: 4
    resources:
        nodes=1,
        mem_mb=16000,
        runtime=240,
    input:
        resources=rules.preparePRSMethodResources.output.ready,
        env=rules.preparePRSInputs.output.env,
    output:
        done=PRS_METHOD_RUN_DIR / "single_prsice.done",
    params:
        enabled=prs_method_enabled("single_prsice"),
        out_dir=PRS_METHOD_RUN_DIR / "single_prsice",
        script=Path(workflow.basedir) / "scripts" / "PRSice.R",
        extra=get_method_extra_args("single_prsice"),
    shell:
        """
        if [[ "{params.enabled}" != "True" ]]; then
            echo "Method single_prsice not enabled, skipping..."
            touch {output.done}
            exit 0
        fi
        source {input.env}
        mkdir -p {params.out_dir}
        Rscript {params.script} \
            --prsice /opt/conda/envs/singlePRS/bin/PRSice \
            --base "$target_sumstats_file" \
            --target "$study_sample_plink" \
            --pheno "$target_study_pheno_file" \
            --out {params.out_dir}/prsice \
            {params.extra} \
            > {log} 2>&1
        touch {output.done}
        """


rule runSingleAncestryPRSCS:
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/prs:latest"
    log:
        OUT_DIR / "logs" / "runSingleAncestryPRSCS.log",
    threads: 4
    resources:
        nodes=1,
        mem_mb=32000,
        runtime=720,
    input:
        resources=rules.preparePRSMethodResources.output.ready,
        env=rules.preparePRSInputs.output.env,
    output:
        done=PRS_METHOD_RUN_DIR / "single_prscs.done",
    params:
        enabled=prs_method_enabled("single_prscs"),
        out_dir=PRS_METHOD_RUN_DIR / "single_prscs",
        extra=get_method_extra_args("single_prscs"),
    shell:
        """
        if [[ "{params.enabled}" != "True" ]]; then
            echo "Method single_prscs not enabled, skipping..."
            touch {output.done}
            exit 0
        fi
        source {input.env}
        mkdir -p {params.out_dir}
        python /usr/local/bin/PRSCS.py \
            --target-sumstats "$target_sumstats_file" \
            --target-pheno "$target_study_pheno_file" \
            --out-dir {params.out_dir} \
            {params.extra} \
            > {log} 2>&1
        touch {output.done}
        """


rule runSingleAncestryLDpred2:
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/prs:latest"
    log:
        OUT_DIR / "logs" / "runSingleAncestryLDpred2.log",
    threads: 4
    resources:
        nodes=1,
        mem_mb=64000,
        runtime=720,
    input:
        resources=rules.preparePRSMethodResources.output.ready,
        env=rules.preparePRSInputs.output.env,
    output:
        done=PRS_METHOD_RUN_DIR / "single_ldpred2.done",
    params:
        enabled=prs_method_enabled("single_ldpred2"),
        out_dir=PRS_METHOD_RUN_DIR / "single_ldpred2",
        script=Path(workflow.basedir) / "scripts" / "ldpred2.R",
        extra=get_method_extra_args("single_ldpred2"),
    shell:
        """
        if [[ "{params.enabled}" != "True" ]]; then
            echo "Method single_ldpred2 not enabled, skipping..."
            touch {output.done}
            exit 0
        fi
        source {input.env}
        mkdir -p {params.out_dir}
        export OMP_NUM_THREADS=1
        Rscript {params.script} \
            --ss "$target_sumstats_file" \
            --anc_bed "$study_sample_plink.bed" \
            --bim "$study_sample_plink.bim" \
            --out {params.out_dir}/ldpred2 \
            {params.extra} \
            > {log} 2>&1
        touch {output.done}
        """


rule runMultiAncestryCTSLEB:
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/prs:latest"
    log:
        OUT_DIR / "logs" / "runMultiAncestryCTSLEB.log",
    threads: 4
    resources:
        nodes=1,
        mem_mb=32000,
        runtime=720,
    input:
        resources=rules.preparePRSMethodResources.output.ready,
        env=rules.preparePRSInputs.output.env,
    output:
        done=PRS_METHOD_RUN_DIR / "multi_ctsleb.done",
    params:
        enabled=prs_method_enabled("multi_ctsleb"),
        out_dir=PRS_METHOD_RUN_DIR / "multi_ctsleb",
        extra=get_method_extra_args("multi_ctsleb"),
    shell:
        """
        if [[ "{params.enabled}" != "True" ]]; then
            echo "Method multi_ctsleb not enabled, skipping..."
            touch {output.done}
            exit 0
        fi
        source {input.env}
        mkdir -p {params.out_dir}
        Rscript /usr/local/bin/ctsleb_multi.R \
            --target-sumstats "$target_sumstats_file" \
            --training-sumstats "$training_sumstats_file" \
            --target-bed "$study_sample_plink.bed" \
            --target-bim "$study_sample_plink.bim" \
            --target-fam "$study_sample_plink.fam" \
            --target-pheno "$target_study_pheno_file" \
            --anc2-bed "$study_sample_plink_anc2.bed" \
            --anc2-bim "$study_sample_plink_anc2.bim" \
            --anc2-fam "$study_sample_plink_anc2.fam" \
            --anc2-pheno "$training_study_pheno_file" \
            --out-dir {params.out_dir} \
            {params.extra} \
            > {log} 2>&1
        touch {output.done}
        """


rule runMultiAncestryPRSCSx:
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/prs:latest"
    log:
        OUT_DIR / "logs" / "runMultiAncestryPRSCSx.log",
    threads: 4
    resources:
        nodes=1,
        mem_mb=32000,
        runtime=720,
    input:
        resources=rules.preparePRSMethodResources.output.ready,
        env=rules.preparePRSInputs.output.env,
    output:
        done=PRS_METHOD_RUN_DIR / "multi_prscsx.done",
    params:
        enabled=prs_method_enabled("multi_prscsx"),
        out_dir=PRS_METHOD_RUN_DIR / "multi_prscsx",
        extra=get_method_extra_args("multi_prscsx"),
    shell:
        """
        if [[ "{params.enabled}" != "True" ]]; then
            echo "Method multi_prscsx not enabled, skipping..."
            touch {output.done}
            exit 0
        fi
        source {input.env}
        mkdir -p {params.out_dir}
        python /usr/local/bin/PRSCSx.py \
            --target-sumstats "$target_sumstats_file" \
            --training-sumstats "$training_sumstats_file" \
            --target-pheno "$target_study_pheno_file" \
            --anc2-pheno "$training_study_pheno_file" \
            --out-dir {params.out_dir} \
            {params.extra} \
            > {log} 2>&1
        touch {output.done}
        """


rule runMultiAncestryLDpred2:
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/prs:latest"
    log:
        OUT_DIR / "logs" / "runMultiAncestryLDpred2.log",
    threads: 4
    resources:
        nodes=1,
        mem_mb=64000,
        runtime=720,
    input:
        resources=rules.preparePRSMethodResources.output.ready,
        env=rules.preparePRSInputs.output.env,
    output:
        done=PRS_METHOD_RUN_DIR / "multi_ldpred2.done",
    params:
        enabled=prs_method_enabled("multi_ldpred2"),
        out_dir=PRS_METHOD_RUN_DIR / "multi_ldpred2",
        extra=get_method_extra_args("multi_ldpred2"),
    shell:
        """
        if [[ "{params.enabled}" != "True" ]]; then
            echo "Method multi_ldpred2 not enabled, skipping..."
            touch {output.done}
            exit 0
        fi
        source {input.env}
        mkdir -p {params.out_dir}
        Rscript /usr/local/bin/ldpred2_multi.R \
            --target-sumstats "$target_sumstats_file" \
            --training-sumstats "$training_sumstats_file" \
            --target-bed "$study_sample_plink.bed" \
            --target-bim "$study_sample_plink.bim" \
            --target-fam "$study_sample_plink.fam" \
            --target-pheno "$target_study_pheno_file" \
            --anc2-bed "$study_sample_plink_anc2.bed" \
            --anc2-bim "$study_sample_plink_anc2.bim" \
            --anc2-fam "$study_sample_plink_anc2.fam" \
            --anc2-pheno "$training_study_pheno_file" \
            --out-dir {params.out_dir} \
            {params.extra} \
            > {log} 2>&1
        touch {output.done}
        """


rule runMultiAncestryPROSPER:
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/prs:latest"
    log:
        OUT_DIR / "logs" / "runMultiAncestryPROSPER.log",
    threads: 4
    resources:
        nodes=1,
        mem_mb=64000,
        runtime=720,
    input:
        resources=rules.preparePRSMethodResources.output.ready,
        env=rules.preparePRSInputs.output.env,
    output:
        done=PRS_METHOD_RUN_DIR / "multi_prosper.done",
    params:
        enabled=prs_method_enabled("multi_prosper"),
        out_dir=PRS_METHOD_RUN_DIR / "multi_prosper",
        extra=get_method_extra_args("multi_prosper"),
    shell:
        """
        if [[ "{params.enabled}" != "True" ]]; then
            echo "Method multi_prosper not enabled, skipping..."
            touch {output.done}
            exit 0
        fi
        source {input.env}
        mkdir -p {params.out_dir}
        python /usr/local/bin/PROSPER.py \
            --target-sumstats "$target_sumstats_file" \
            --training-sumstats "$training_sumstats_file" \
            --target-pheno "$target_study_pheno_file" \
            --anc2-pheno "$training_study_pheno_file" \
            --out-dir {params.out_dir} \
            {params.extra} \
            > {log} 2>&1
        touch {output.done}
        """


rule runMultiAncestrySDPRS:
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/prs:latest"
    log:
        OUT_DIR / "logs" / "runMultiAncestrySDPRS.log",
    threads: 4
    resources:
        nodes=1,
        mem_mb=64000,
        runtime=720,
    input:
        resources=rules.preparePRSMethodResources.output.ready,
        env=rules.preparePRSInputs.output.env,
    output:
        done=PRS_METHOD_RUN_DIR / "multi_sdprs.done",
    params:
        enabled=prs_method_enabled("multi_sdprs"),
        out_dir=PRS_METHOD_RUN_DIR / "multi_sdprs",
        extra=get_method_extra_args("multi_sdprs"),
    shell:
        """
        if [[ "{params.enabled}" != "True" ]]; then
            echo "Method multi_sdprs not enabled, skipping..."
            touch {output.done}
            exit 0
        fi
        source {input.env}
        mkdir -p {params.out_dir}
        python /usr/local/bin/SDPRX.py \
            --target-sumstats "$target_sumstats_file" \
            --training-sumstats "$training_sumstats_file" \
            --target-pheno "$target_study_pheno_file" \
            --anc2-pheno "$training_study_pheno_file" \
            --out-dir {params.out_dir} \
            {params.extra} \
            > {log} 2>&1
        touch {output.done}
        """


rule runAllEnabledPRS:
    """Run all PRS methods that are enabled in config."""
    input:
        resources=rules.preparePRSMethodResources.output.ready,
        env=rules.preparePRSInputs.output.env,
        single_ct_done=PRS_METHOD_RUN_DIR / "single_ct.done" if prs_method_enabled("single_ct") else [],
        single_prsice_done=PRS_METHOD_RUN_DIR / "single_prsice.done" if prs_method_enabled("single_prsice") else [],
        single_prscs_done=PRS_METHOD_RUN_DIR / "single_prscs.done" if prs_method_enabled("single_prscs") else [],
        single_ldpred2_done=PRS_METHOD_RUN_DIR / "single_ldpred2.done" if prs_method_enabled("single_ldpred2") else [],
        single_lassosum2_done=PRS_METHOD_RUN_DIR / "single_lassosum2.done" if prs_method_enabled("single_lassosum2") else [],
        multi_ctsleb_done=PRS_METHOD_RUN_DIR / "multi_ctsleb.done" if prs_method_enabled("multi_ctsleb") else [],
        multi_prscsx_done=PRS_METHOD_RUN_DIR / "multi_prscsx.done" if prs_method_enabled("multi_prscsx") else [],
        multi_ldpred2_done=PRS_METHOD_RUN_DIR / "multi_ldpred2.done" if prs_method_enabled("multi_ldpred2") else [],
        multi_prosper_done=PRS_METHOD_RUN_DIR / "multi_prosper.done" if prs_method_enabled("multi_prosper") else [],
        multi_sdprs_done=PRS_METHOD_RUN_DIR / "multi_sdprs.done" if prs_method_enabled("multi_sdprs") else [],
    output:
        done=OUT_DIR / "prs_all_completed.done",
    log:
        OUT_DIR / "logs" / "runAllEnabledPRS.log",
    params:
        enabled_methods=get_enabled_prs_methods(),
    shell:
        """
        echo "All enabled PRS methods completed." > {log} 2>&1
        echo "Enabled methods: {params.enabled_methods}" >> {log} 2>&1
        touch {output.done}
        """
