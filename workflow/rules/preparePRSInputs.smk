PRS_CONFIG = config.get("prsPipeline", {})
PRS_SIM_CONFIG = config.get("phenotypeSimulation", {})

# Absorbed prs_pipeline method engine (see workflow/scripts/prs_pipeline/VENDORED.md).
# Default = the vendored in-repo copy; prs_pipeline_path may override it (e.g. for
# legacy MSI configs). prs_pipeline_ref records the upstream source commit (pin).
PRS_PIPELINE_PATH = Path(
    config.get("prs_pipeline_path")
    or str(Path(workflow.basedir) / "scripts" / "prs_pipeline")
)
PRS_PIPELINE_PIN = config.get(
    "prs_pipeline_ref",
    "f71cf4f1031ff6231bc9f8abc9f96ebdc6f17bdf",  # sandbox_multi_pheno (2026-08-20)
)

PRS_ANC1 = PRS_SIM_CONFIG.get("ancestries", ["AFR", "EUR"])[0]
PRS_ANC2 = PRS_SIM_CONFIG.get("ancestries", ["AFR", "EUR"])[1]
PRS_SIM_DIR = OUT_DIR / "simulations" / f"{PRS_ANC1}_{PRS_ANC2}"
PRS_OUT_DIR = Path(
    PRS_CONFIG.get(
        "generated_input_dir",
        str(OUT_DIR / "prs_inputs" / f"{PRS_ANC1}_{PRS_ANC2}"),
    )
)


rule checkPRSPipelinePath:
    log:
        OUT_DIR / "logs" / "checkPRSPipelinePath.log",
    threads: 1
    resources:
        nodes=1,
        mem_mb=2000,
        runtime=30,
    output:
        checked=OUT_DIR / ".prs_pipeline_path.checked",
    params:
        repo=PRS_PIPELINE_PATH,
        pin=PRS_PIPELINE_PIN,
    shell:
        """
        set -euo pipefail

        REPO="{params.repo}"
        PIN="{params.pin}"

        if [[ ! -d "$REPO" ]]; then
            echo "ERROR: prs_pipeline tree not found at $REPO (check prs_pipeline_path)" >&2
            exit 1
        fi
        for f in run_single_ancestry_PRS_pipeline.sh \\
                 src/prepare_sumstats.R \\
                 src/run_CT.sh \\
                 src/run_LDpred2.R \\
                 src/run_lassosum2.R \\
                 src/run_PRSice2.sh; do
            if [[ ! -f "$REPO/$f" ]]; then
                echo "ERROR: missing vendored $REPO/$f" >&2
                exit 1
            fi
        done

        echo "prs_pipeline_path=$REPO" > {output.checked}
        echo "prs_pipeline_ref=$PIN" >> {output.checked}
        echo "prs_pipeline_vendored=true" >> {output.checked}
        echo "checked_at=$(date -Iseconds)" >> {output.checked}

        echo "OK: prs_pipeline tree found at $REPO (source pin $PIN)" > "{log}"
        """


rule preparePRSInputs:
    log:
        OUT_DIR / "logs" / "preparePRSInputs.log",
    container:
        "oras://ghcr.io/mainsqu33ze/gdcgenomicsqc/prs:latest"
    conda:
        "../../envs/prs.yml"
    threads: 4
    resources:
        nodes=1,
        mem_mb=16000,
        runtime=120,
    input:
        anc1_bed=PRS_SIM_DIR / f"{PRS_ANC1}_simulation.bed",
        anc1_bim=PRS_SIM_DIR / f"{PRS_ANC1}_simulation.bim",
        anc1_fam=PRS_SIM_DIR / f"{PRS_ANC1}_simulation.fam",
        anc2_bed=PRS_SIM_DIR / f"{PRS_ANC2}_simulation.bed",
        anc2_bim=PRS_SIM_DIR / f"{PRS_ANC2}_simulation.bim",
        anc2_fam=PRS_SIM_DIR / f"{PRS_ANC2}_simulation.fam",
    output:
        target_sumstats=PRS_OUT_DIR / "gwas" / "target_sumstats.txt",
        training_sumstats=PRS_OUT_DIR / "gwas" / "training_sumstats.txt",
        target_single_sumstats=PRS_OUT_DIR / "gwas" / "target_sumstats_singlePRS.txt",
        training_single_sumstats=PRS_OUT_DIR / "gwas" / "training_sumstats_singlePRS.txt",
        target_gwas_pheno=PRS_OUT_DIR / "metadata" / f"{PRS_ANC1}_gwas.pheno",
        target_study_pheno=PRS_OUT_DIR / "metadata" / f"{PRS_ANC1}_study.pheno",
        training_gwas_pheno=PRS_OUT_DIR / "metadata" / f"{PRS_ANC2}_gwas.pheno",
        training_study_pheno=PRS_OUT_DIR / "metadata" / f"{PRS_ANC2}_study.pheno",
        study_bed=PRS_OUT_DIR / "anc1_plink_files" / f"{PRS_ANC1}_simulation_study_sample.bed",
        study_bim=PRS_OUT_DIR / "anc1_plink_files" / f"{PRS_ANC1}_simulation_study_sample.bim",
        study_fam=PRS_OUT_DIR / "anc1_plink_files" / f"{PRS_ANC1}_simulation_study_sample.fam",
        study_anc2_bed=PRS_OUT_DIR / "anc2_plink_files" / f"{PRS_ANC2}_simulation_study_sample.bed",
        study_anc2_bim=PRS_OUT_DIR / "anc2_plink_files" / f"{PRS_ANC2}_simulation_study_sample.bim",
        study_anc2_fam=PRS_OUT_DIR / "anc2_plink_files" / f"{PRS_ANC2}_simulation_study_sample.fam",
        env=PRS_OUT_DIR / "prs_inputs.env",
        prscsx_config=PRS_OUT_DIR / "prs_prscsx_generated.conf",
        single_config=PRS_OUT_DIR / f"prs_single_ancestry_{PRS_ANC1}_generated.conf",
    params:
        sim_dir=PRS_SIM_DIR,
        out_dir=PRS_OUT_DIR,
        anc1=PRS_ANC1,
        anc2=PRS_ANC2,
        phenotype_index=PRS_CONFIG.get("phenotype_index", 1),
        gwas_fraction=PRS_CONFIG.get("gwas_fraction", 0.5),
        seed=PRS_CONFIG.get("seed", 42),
        plink2=PRS_CONFIG.get("path_plink2", "plink2"),
        prs_pipeline_dir=PRS_PIPELINE_PATH,
        script=Path(workflow.basedir) / "scripts" / "prepare_prs_inputs.sh",
    shell:
        """
        bash {params.script} \
            --sim-dir {params.sim_dir} \
            --out-dir {params.out_dir} \
            --anc1 {params.anc1} \
            --anc2 {params.anc2} \
            --phenotype-index {params.phenotype_index} \
            --gwas-fraction {params.gwas_fraction} \
            --seed {params.seed} \
            --plink2-bin {params.plink2} \
            --prs-pipeline-dir {params.prs_pipeline_dir} \
            > {log} 2>&1
        """


rule runSingleAncestryPRS:
    log:
        OUT_DIR / "logs" / f"runSingleAncestryPRS_{PRS_ANC1}.log",
    container:
        "oras://ghcr.io/mainsqu33ze/gdcgenomicsqc/prs:latest"
    conda:
        "../../envs/prs.yml"
    threads: 4
    resources:
        nodes=1,
        mem_mb=16000,
        runtime=240,
    input:
        check=rules.checkPRSPipelinePath.output.checked,
        config=rules.preparePRSInputs.output.single_config,
        target_sumstats=rules.preparePRSInputs.output.target_single_sumstats,
        study_bed=rules.preparePRSInputs.output.study_bed,
        study_bim=rules.preparePRSInputs.output.study_bim,
        study_fam=rules.preparePRSInputs.output.study_fam,
    output:
        done=PRS_OUT_DIR / f"single_ancestry_{PRS_ANC1}.done",
    params:
        script=PRS_CONFIG.get(
            "single_ancestry_script",
            str(PRS_PIPELINE_PATH / "run_single_ancestry_PRS_pipeline.sh"),
        ),
        flags=PRS_CONFIG.get("single_ancestry_flags", "-c -l -s -P"),
    shell:
        """
        set -euo pipefail

        echo "Running single-ancestry PRS pipeline" > {log}
        echo "Script: {params.script}" >> {log}
        echo "Config: {input.config}" >> {log}
        echo "Flags: {params.flags}" >> {log}

        if [[ ! -f "{params.script}" ]]; then
            echo "Missing single-ancestry PRS script: {params.script}" >> {log}
            exit 1
        fi

        # PRSice2: drive the vendored PRSice.R + static PRSice_linux directly (same
        # pattern as runSingleAncestryPRSice; matches both the conda-PRSice `prs`
        # image and any image lacking a `PRSice` command).
        export PRSICE_CMD="Rscript {PRS_PIPELINE_PATH}/src/PRSice.R --prsice {PRS_PIPELINE_PATH}/src/PRSice_linux"

        bash {params.script} {params.flags} -C {input.config} >> {log} 2>&1

        touch {output.done}
        """
