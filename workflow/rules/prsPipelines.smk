import shlex

PRS_METHODS_CONFIG = config.get("prsMethods", {})
PRS_RESOURCE_DIR = Path(
    PRS_METHODS_CONFIG.get(
        "resource_dir",
        str(OUT_DIR.parent / "prs_resources"),
    )
)
PRS_METHOD_RUN_DIR = PRS_OUT_DIR / "method_runs"
PRS_DOWNLOAD_SOFTWARE_FLAG = (
    "--download-software" if PRS_METHODS_CONFIG.get("download_software", False) else ""
)


def prs_method_command(method):
    return PRS_METHODS_CONFIG.get(method, {}).get("command", "")


def prs_method_command_quoted(method):
    return shlex.quote(prs_method_command(method))


def prs_method_extra_args(method):
    method_config = PRS_METHODS_CONFIG.get(method, {})
    args = []
    for key, flag in (
        ("ld_ref_dir", "--ld-ref-dir"),
        ("ld_ref_prefix", "--ld-ref-prefix"),
        ("ld_matrix_dir", "--ld-matrix-dir"),
        ("software_dir", "--software-dir"),
    ):
        value = method_config.get(key)
        if value:
            args.extend([flag, shlex.quote(str(value))])
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
    params:
        resource_dir=PRS_RESOURCE_DIR,
        script=Path(workflow.basedir) / "scripts" / "download_prs_resources.sh",
        prscsx_ref=PRS_CONFIG.get(
            "path_ref_dir",
            PRS_METHODS_CONFIG.get("multi_prscsx", {}).get("ld_ref_dir", ""),
        ),
        plink2=PRS_CONFIG.get("path_plink2", PRS_METHODS_CONFIG.get("plink2", "")),
        download_software=PRS_DOWNLOAD_SOFTWARE_FLAG,
    shell:
        """
        bash {params.script} \
            --resource-dir {params.resource_dir} \
            --prscsx-ref-dir {params.prscsx_ref} \
            --plink2 {params.plink2} \
            {params.download_software} \
            > {log} 2>&1
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
        sumstats=rules.preparePRSInputs.output.target_sumstats,
        bed=rules.preparePRSInputs.output.study_bed,
        bim=rules.preparePRSInputs.output.study_bim,
        fam=rules.preparePRSInputs.output.study_fam,
        pheno=rules.preparePRSInputs.output.target_study_pheno,
    output:
        done=PRS_METHOD_RUN_DIR / "single_ct.done",
    params:
        method="single_ct",
        command=prs_method_command_quoted("single_ct"),
        extra=prs_method_extra_args("single_ct"),
        out_dir=PRS_METHOD_RUN_DIR / "single_ct",
        script=Path(workflow.basedir) / "scripts" / "run_prs_pipeline_adapter.sh",
    shell:
        """
        PRS_METHOD_COMMAND={params.command} bash {params.script} \
            --method {params.method} \
            --prs-inputs-env {input.env} \
            --resource-dir {PRS_RESOURCE_DIR} \
            --out-dir {params.out_dir} \
            {params.extra} \
            --done {output.done} \
            > {log} 2>&1
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
        sumstats=rules.preparePRSInputs.output.target_sumstats,
        bed=rules.preparePRSInputs.output.study_bed,
        bim=rules.preparePRSInputs.output.study_bim,
        fam=rules.preparePRSInputs.output.study_fam,
        pheno=rules.preparePRSInputs.output.target_study_pheno,
    output:
        done=PRS_METHOD_RUN_DIR / "single_prsice.done",
    params:
        method="single_prsice",
        command=prs_method_command_quoted("single_prsice"),
        extra=prs_method_extra_args("single_prsice"),
        out_dir=PRS_METHOD_RUN_DIR / "single_prsice",
        script=Path(workflow.basedir) / "scripts" / "run_prs_pipeline_adapter.sh",
    shell:
        """
        PRS_METHOD_COMMAND={params.command} bash {params.script} \
            --method {params.method} \
            --prs-inputs-env {input.env} \
            --resource-dir {PRS_RESOURCE_DIR} \
            --out-dir {params.out_dir} \
            {params.extra} \
            --done {output.done} \
            > {log} 2>&1
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
        sumstats=rules.preparePRSInputs.output.target_sumstats,
        bim=rules.preparePRSInputs.output.study_bim,
        pheno=rules.preparePRSInputs.output.target_study_pheno,
    output:
        done=PRS_METHOD_RUN_DIR / "single_prscs.done",
    params:
        method="single_prscs",
        command=prs_method_command_quoted("single_prscs"),
        extra=prs_method_extra_args("single_prscs"),
        out_dir=PRS_METHOD_RUN_DIR / "single_prscs",
        script=Path(workflow.basedir) / "scripts" / "run_prs_pipeline_adapter.sh",
    shell:
        """
        PRS_METHOD_COMMAND={params.command} bash {params.script} \
            --method {params.method} \
            --prs-inputs-env {input.env} \
            --resource-dir {PRS_RESOURCE_DIR} \
            --out-dir {params.out_dir} \
            {params.extra} \
            --done {output.done} \
            > {log} 2>&1
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
        sumstats=rules.preparePRSInputs.output.target_sumstats,
        bed=rules.preparePRSInputs.output.study_bed,
        bim=rules.preparePRSInputs.output.study_bim,
        fam=rules.preparePRSInputs.output.study_fam,
        pheno=rules.preparePRSInputs.output.target_study_pheno,
    output:
        done=PRS_METHOD_RUN_DIR / "single_ldpred2.done",
    params:
        method="single_ldpred2",
        command=prs_method_command_quoted("single_ldpred2"),
        extra=prs_method_extra_args("single_ldpred2"),
        out_dir=PRS_METHOD_RUN_DIR / "single_ldpred2",
        script=Path(workflow.basedir) / "scripts" / "run_prs_pipeline_adapter.sh",
    shell:
        """
        PRS_METHOD_COMMAND={params.command} bash {params.script} \
            --method {params.method} \
            --prs-inputs-env {input.env} \
            --resource-dir {PRS_RESOURCE_DIR} \
            --out-dir {params.out_dir} \
            {params.extra} \
            --done {output.done} \
            > {log} 2>&1
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
        sumstats=rules.preparePRSInputs.output.target_sumstats,
        bed=rules.preparePRSInputs.output.study_bed,
        bim=rules.preparePRSInputs.output.study_bim,
        fam=rules.preparePRSInputs.output.study_fam,
        pheno=rules.preparePRSInputs.output.target_study_pheno,
    output:
        done=PRS_METHOD_RUN_DIR / "single_lassosum2.done",
    params:
        method="single_lassosum2",
        command=prs_method_command_quoted("single_lassosum2"),
        extra=prs_method_extra_args("single_lassosum2"),
        out_dir=PRS_METHOD_RUN_DIR / "single_lassosum2",
        script=Path(workflow.basedir) / "scripts" / "run_prs_pipeline_adapter.sh",
    shell:
        """
        PRS_METHOD_COMMAND={params.command} bash {params.script} \
            --method {params.method} \
            --prs-inputs-env {input.env} \
            --resource-dir {PRS_RESOURCE_DIR} \
            --out-dir {params.out_dir} \
            {params.extra} \
            --done {output.done} \
            > {log} 2>&1
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
        target_sumstats=rules.preparePRSInputs.output.target_sumstats,
        training_sumstats=rules.preparePRSInputs.output.training_sumstats,
        study_bed=rules.preparePRSInputs.output.study_bed,
        study_anc2_bed=rules.preparePRSInputs.output.study_anc2_bed,
        study_pheno=rules.preparePRSInputs.output.target_study_pheno,
        study_anc2_pheno=rules.preparePRSInputs.output.training_study_pheno,
    output:
        done=PRS_METHOD_RUN_DIR / "multi_ctsleb.done",
    params:
        method="multi_ctsleb",
        command=prs_method_command_quoted("multi_ctsleb"),
        extra=prs_method_extra_args("multi_ctsleb"),
        out_dir=PRS_METHOD_RUN_DIR / "multi_ctsleb",
        script=Path(workflow.basedir) / "scripts" / "run_prs_pipeline_adapter.sh",
    shell:
        """
        PRS_METHOD_COMMAND={params.command} bash {params.script} \
            --method {params.method} \
            --prs-inputs-env {input.env} \
            --resource-dir {PRS_RESOURCE_DIR} \
            --out-dir {params.out_dir} \
            {params.extra} \
            --done {output.done} \
            > {log} 2>&1
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
        config=rules.preparePRSInputs.output.prscsx_config,
        target_sumstats=rules.preparePRSInputs.output.target_sumstats,
        training_sumstats=rules.preparePRSInputs.output.training_sumstats,
        study_pheno=rules.preparePRSInputs.output.target_study_pheno,
        study_anc2_pheno=rules.preparePRSInputs.output.training_study_pheno,
    output:
        done=PRS_METHOD_RUN_DIR / "multi_prscsx.done",
    params:
        method="multi_prscsx",
        command=prs_method_command_quoted("multi_prscsx"),
        extra=prs_method_extra_args("multi_prscsx"),
        out_dir=PRS_METHOD_RUN_DIR / "multi_prscsx",
        script=Path(workflow.basedir) / "scripts" / "run_prs_pipeline_adapter.sh",
    shell:
        """
        PRS_METHOD_COMMAND={params.command} bash {params.script} \
            --method {params.method} \
            --prs-inputs-env {input.env} \
            --resource-dir {PRS_RESOURCE_DIR} \
            --out-dir {params.out_dir} \
            {params.extra} \
            --done {output.done} \
            > {log} 2>&1
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
        target_sumstats=rules.preparePRSInputs.output.target_sumstats,
        training_sumstats=rules.preparePRSInputs.output.training_sumstats,
        study_bed=rules.preparePRSInputs.output.study_bed,
        study_anc2_bed=rules.preparePRSInputs.output.study_anc2_bed,
        study_pheno=rules.preparePRSInputs.output.target_study_pheno,
        study_anc2_pheno=rules.preparePRSInputs.output.training_study_pheno,
    output:
        done=PRS_METHOD_RUN_DIR / "multi_ldpred2.done",
    params:
        method="multi_ldpred2",
        command=prs_method_command_quoted("multi_ldpred2"),
        extra=prs_method_extra_args("multi_ldpred2"),
        out_dir=PRS_METHOD_RUN_DIR / "multi_ldpred2",
        script=Path(workflow.basedir) / "scripts" / "run_prs_pipeline_adapter.sh",
    shell:
        """
        PRS_METHOD_COMMAND={params.command} bash {params.script} \
            --method {params.method} \
            --prs-inputs-env {input.env} \
            --resource-dir {PRS_RESOURCE_DIR} \
            --out-dir {params.out_dir} \
            {params.extra} \
            --done {output.done} \
            > {log} 2>&1
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
        target_sumstats=rules.preparePRSInputs.output.target_sumstats,
        training_sumstats=rules.preparePRSInputs.output.training_sumstats,
        study_pheno=rules.preparePRSInputs.output.target_study_pheno,
        study_anc2_pheno=rules.preparePRSInputs.output.training_study_pheno,
    output:
        done=PRS_METHOD_RUN_DIR / "multi_prosper.done",
    params:
        method="multi_prosper",
        command=prs_method_command_quoted("multi_prosper"),
        extra=prs_method_extra_args("multi_prosper"),
        out_dir=PRS_METHOD_RUN_DIR / "multi_prosper",
        script=Path(workflow.basedir) / "scripts" / "run_prs_pipeline_adapter.sh",
    shell:
        """
        PRS_METHOD_COMMAND={params.command} bash {params.script} \
            --method {params.method} \
            --prs-inputs-env {input.env} \
            --resource-dir {PRS_RESOURCE_DIR} \
            --out-dir {params.out_dir} \
            {params.extra} \
            --done {output.done} \
            > {log} 2>&1
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
        target_sumstats=rules.preparePRSInputs.output.target_sumstats,
        training_sumstats=rules.preparePRSInputs.output.training_sumstats,
        study_pheno=rules.preparePRSInputs.output.target_study_pheno,
        study_anc2_pheno=rules.preparePRSInputs.output.training_study_pheno,
    output:
        done=PRS_METHOD_RUN_DIR / "multi_sdprs.done",
    params:
        method="multi_sdprs",
        command=prs_method_command_quoted("multi_sdprs"),
        extra=prs_method_extra_args("multi_sdprs"),
        out_dir=PRS_METHOD_RUN_DIR / "multi_sdprs",
        script=Path(workflow.basedir) / "scripts" / "run_prs_pipeline_adapter.sh",
    shell:
        """
        PRS_METHOD_COMMAND={params.command} bash {params.script} \
            --method {params.method} \
            --prs-inputs-env {input.env} \
            --resource-dir {PRS_RESOURCE_DIR} \
            --out-dir {params.out_dir} \
            {params.extra} \
            --done {output.done} \
            > {log} 2>&1
        """
