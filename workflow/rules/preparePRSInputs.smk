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
        target_tuning_keep=PRS_OUT_DIR / "metadata" / f"{PRS_ANC1}_tuning.keep",
        target_testing_keep=PRS_OUT_DIR / "metadata" / f"{PRS_ANC1}_testing.keep",
        target_tuning_pheno=PRS_OUT_DIR / "metadata" / f"{PRS_ANC1}_tuning.pheno",
        target_testing_pheno=PRS_OUT_DIR / "metadata" / f"{PRS_ANC1}_testing.pheno",
        training_gwas_pheno=PRS_OUT_DIR / "metadata" / f"{PRS_ANC2}_gwas.pheno",
        training_study_pheno=PRS_OUT_DIR / "metadata" / f"{PRS_ANC2}_study.pheno",
        study_bed=PRS_OUT_DIR / "anc1_plink_files" / f"{PRS_ANC1}_simulation_study_sample.bed",
        study_bim=PRS_OUT_DIR / "anc1_plink_files" / f"{PRS_ANC1}_simulation_study_sample.bim",
        study_fam=PRS_OUT_DIR / "anc1_plink_files" / f"{PRS_ANC1}_simulation_study_sample.fam",
        target_ld_reference_bed=PRS_OUT_DIR / "anc1_plink_files" / f"{PRS_ANC1}_simulation_ld_reference.bed",
        target_ld_reference_bim=PRS_OUT_DIR / "anc1_plink_files" / f"{PRS_ANC1}_simulation_ld_reference.bim",
        target_ld_reference_fam=PRS_OUT_DIR / "anc1_plink_files" / f"{PRS_ANC1}_simulation_ld_reference.fam",
        ld_study_bed=PRS_OUT_DIR / "anc1_plink_files" / f"{PRS_ANC1}_simulation_study_sample_ldreduced.bed",
        ld_study_bim=PRS_OUT_DIR / "anc1_plink_files" / f"{PRS_ANC1}_simulation_study_sample_ldreduced.bim",
        ld_study_fam=PRS_OUT_DIR / "anc1_plink_files" / f"{PRS_ANC1}_simulation_study_sample_ldreduced.fam",
        ld_variant_ids=PRS_OUT_DIR / "metadata" / f"{PRS_ANC1}_ld_variant_ids.txt",
        study_anc2_bed=PRS_OUT_DIR / "anc2_plink_files" / f"{PRS_ANC2}_simulation_study_sample.bed",
        study_anc2_bim=PRS_OUT_DIR / "anc2_plink_files" / f"{PRS_ANC2}_simulation_study_sample.bim",
        study_anc2_fam=PRS_OUT_DIR / "anc2_plink_files" / f"{PRS_ANC2}_simulation_study_sample.fam",
        aux_ld_reference_bed=PRS_OUT_DIR / "anc2_plink_files" / f"{PRS_ANC2}_simulation_ld_reference.bed",
        aux_ld_reference_bim=PRS_OUT_DIR / "anc2_plink_files" / f"{PRS_ANC2}_simulation_ld_reference.bim",
        aux_ld_reference_fam=PRS_OUT_DIR / "anc2_plink_files" / f"{PRS_ANC2}_simulation_ld_reference.fam",
        target_analysis_variant_extract=PRS_OUT_DIR / "metadata" / f"{PRS_ANC1}_analysis_variants.extract",
        aux_analysis_variant_extract=PRS_OUT_DIR / "metadata" / f"{PRS_ANC2}_analysis_variants.extract",
        cohort_allocation=PRS_OUT_DIR / "metadata" / "cohort_allocation.tsv",
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
        gwas_n=PRS_CONFIG.get("gwas_n", 0),
        ld_reference_n=PRS_CONFIG.get("ld_reference_n", 0),
        tuning_n=PRS_CONFIG.get("tuning_n", 0),
        testing_n=PRS_CONFIG.get("testing_n", 0),
        analysis_variant_mode=PRS_CONFIG.get("analysis_variant_mode", "none"),
        prscsx_snpinfo=PRS_CONFIG.get("prscsx_snpinfo", ""),
        seed=PRS_CONFIG.get("seed", 42),
        plink1=PRS_CONFIG.get("path_plink1", "plink"),
        plink2=PRS_CONFIG.get("path_plink2", "plink2"),
        prs_pipeline=PRS_CONFIG.get(
            "path_prs_pipeline", str(Path(workflow.basedir) / "scripts" / "prs_pipeline")
        ),
        ld_variant_extract=PRS_CONFIG.get("ld_variant_extract", ""),
        ld_prune_window_kb=PRS_CONFIG.get("ld_prune_window_kb", 500),
        ld_prune_step=PRS_CONFIG.get("ld_prune_step", 50),
        ld_prune_r2=PRS_CONFIG.get("ld_prune_r2", 0.1),
        ld_prune_max_variants=PRS_CONFIG.get("ld_prune_max_variants", 100000),
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
            --gwas-n {params.gwas_n} \
            --ld-reference-n {params.ld_reference_n} \
            --tuning-n {params.tuning_n} \
            --testing-n {params.testing_n} \
            --analysis-variant-mode {params.analysis_variant_mode} \
            --prscsx-snpinfo "{params.prscsx_snpinfo}" \
            --seed {params.seed} \
            --plink1-bin {params.plink1} \
            --plink2-bin {params.plink2} \
            --prs-pipeline-dir {params.prs_pipeline} \
            --ld-variant-extract "{params.ld_variant_extract}" \
            --ld-prune-window-kb {params.ld_prune_window_kb} \
            --ld-prune-step {params.ld_prune_step} \
            --ld-prune-r2 {params.ld_prune_r2} \
            --ld-prune-max-variants {params.ld_prune_max_variants} \
            > {log} 2>&1
        """


rule runSingleAncestryPRS:
    log:
        OUT_DIR / "logs" / f"runSingleAncestryPRS_{PRS_ANC1}.log",
    threads: int(PRS_CONFIG.get("single_ancestry_threads", 4))
    resources:
        nodes=1,
        mem_mb=int(PRS_CONFIG.get("single_ancestry_mem_mb", 64000)),
        runtime=int(PRS_CONFIG.get("single_ancestry_runtime", 720)),
    input:
        config=rules.preparePRSInputs.output.single_config,
        target_sumstats=rules.preparePRSInputs.output.target_single_sumstats,
        study_bed=rules.preparePRSInputs.output.study_bed,
        study_bim=rules.preparePRSInputs.output.study_bim,
        study_fam=rules.preparePRSInputs.output.study_fam,
        ld_study_bed=rules.preparePRSInputs.output.ld_study_bed,
        ld_study_bim=rules.preparePRSInputs.output.ld_study_bim,
        ld_study_fam=rules.preparePRSInputs.output.ld_study_fam,
        tuning_keep=rules.preparePRSInputs.output.target_tuning_keep,
        testing_keep=rules.preparePRSInputs.output.target_testing_keep,
        tuning_pheno=rules.preparePRSInputs.output.target_tuning_pheno,
        testing_pheno=rules.preparePRSInputs.output.target_testing_pheno,
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
