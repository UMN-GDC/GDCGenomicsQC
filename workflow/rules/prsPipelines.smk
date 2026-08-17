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

CTSLEB_INPUT_DIR = Path(config.get("prsPipeline", {}).get("generated_input_dir", ""))
CTSLEB_PHENO_PREFIX = config.get("prsMethods", {}).get("multi_ctsleb", {}).get(
    "pheno_prefix",
    f"{PRS_ANC1}_HM3",
)

rule prepareCTSLEBSumstats:
    log:
        OUT_DIR / "logs" / "prepareCTSLEBSumstats.log",
    input:
        target_ss=rules.preparePRSInputs.output.target_sumstats,
        training_ss=rules.preparePRSInputs.output.training_sumstats,
        target_bim=rules.preparePRSInputs.output.study_bim,
        training_bim=rules.preparePRSInputs.output.study_anc2_bim,
    output:
        target_ss=CTSLEB_INPUT_DIR / "gwas" / "target_sumstats_CTSLEB.txt",
        training_ss=CTSLEB_INPUT_DIR / "gwas" / "training_sumstats_CTSLEB.txt",
    run:
        import pandas as pd

        def make_ctsleb_sumstats(ss_file, bim_file, out_file):
            ss = pd.read_csv(str(ss_file), sep="\t")
            bim = pd.read_csv(
                str(bim_file),
                sep=r"\s+",
                header=None,
                names=["CHR", "SNP", "CM", "BP", "A1_bim", "A2_bim"],
            )

            if "SNP" not in ss.columns:
                raise ValueError(f"Missing SNP column in {ss_file}")

            merged = ss.merge(bim[["SNP", "BP"]], on="SNP", how="left")
            missing_bp = merged["BP"].isna().sum()
            if missing_bp:
                raise ValueError(f"Missing BP for {missing_bp} rows in {out_file}")

            merged["BP"] = merged["BP"].astype(int)
            merged["rs_id"] = merged["SNP"]

            keep = ["CHR", "SNP", "BP", "A1", "A2", "BETA", "SE", "P", "N", "rs_id"]
            missing_cols = [col for col in keep if col not in merged.columns]
            if missing_cols:
                raise ValueError(
                    f"Missing columns in merged file {out_file}: {', '.join(missing_cols)}"
                )

            merged.loc[:, keep].to_csv(str(out_file), sep="\t", index=False)

        make_ctsleb_sumstats(input.target_ss, input.target_bim, output.target_ss)
        make_ctsleb_sumstats(input.training_ss, input.training_bim, output.training_ss)

        with open(str(log), "w", encoding="utf-8") as handle:
            handle.write(f"Wrote: {output.target_ss}\n")
            handle.write(f"Wrote: {output.training_ss}\n")

rule prepareCTSLEBPhenotypes:
    log:
        OUT_DIR / "logs" / "prepareCTSLEBPhenotypes.log",
    threads: 1
    resources:
        nodes=1,
        mem_mb=4000,
        runtime=60,
    input:
        fam=rules.preparePRSInputs.output.study_fam,
    output:
        tuning=CTSLEB_INPUT_DIR / "metadata" / f"{CTSLEB_PHENO_PREFIX}_tuning.pheno",
        validation=CTSLEB_INPUT_DIR / "metadata" / f"{CTSLEB_PHENO_PREFIX}_validation.pheno",
    shell:
        """
        set -euo pipefail

        mkdir -p $(dirname {output.tuning})

        N=$(wc -l < {input.fam})
        N_TUNE=$((N / 2))
        awk -v n_tune="$N_TUNE" -v tune="{output.tuning}" -v valid="{output.validation}" '
            NR <= n_tune {{ print $6 > tune; next }}
            {{ print $6 > valid }}
        ' {input.fam}

        echo "Created CTSLEB phenotype splits:" > {log}
        wc -l {output.tuning} {output.validation} >> {log}
        """

rule runMultiAncestryCTSLEB:
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
        target_sumstats=rules.prepareCTSLEBSumstats.output.target_ss,
        training_sumstats=rules.prepareCTSLEBSumstats.output.training_ss,
        study_bed=rules.preparePRSInputs.output.study_bed,
        study_anc2_bed=rules.preparePRSInputs.output.study_anc2_bed,
        study_pheno=rules.preparePRSInputs.output.target_study_pheno,
        study_anc2_pheno=rules.preparePRSInputs.output.training_study_pheno,
        tuning_pheno=rules.prepareCTSLEBPhenotypes.output.tuning,
        validation_pheno=rules.prepareCTSLEBPhenotypes.output.validation,
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



rule runMultiAncestrySDPRS:
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
