PRS_CONFIG = config.get("prsPipeline", {})
PRS_SIM_CONFIG = config.get("phenotypeSimulation", {})

PRS_ANC1 = PRS_SIM_CONFIG.get("ancestries", ["AFR", "EUR"])[0]
PRS_ANC2 = PRS_SIM_CONFIG.get("ancestries", ["AFR", "EUR"])[1]
PRS_SIM_DIR = OUT_DIR / "simulations" / f"{PRS_ANC1}_{PRS_ANC2}"
PRS_OUT_DIR = Path(
    PRS_CONFIG.get(
        "generated_input_dir",
        str(OUT_DIR / "prs_inputs" / f"{PRS_ANC1}_{PRS_ANC2}"),
    )
)


rule preparePRSInputs:
    log:
        OUT_DIR / "logs" / "preparePRSInputs.log",
    conda:
        "../../envs/phenotypeSim.yml"
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/prs:latest"
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
            > {log} 2>&1
        """


rule runSingleAncestryPRS:
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/prs:latest"
    log:
        OUT_DIR / "logs" / f"runSingleAncestryPRS_{PRS_ANC1}.log",
    threads: 4
    resources:
        nodes=1,
        mem_mb=16000,
        runtime=240,
    input:
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
            "/projects/standard/gdc/public/prs_methods/scripts/prs_pipeline/run_single_ancestry_PRS_pipeline.sh",
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

        bash {params.script} {params.flags} -C {input.config} >> {log} 2>&1

        touch {output.done}
        """
