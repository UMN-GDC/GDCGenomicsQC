
checkpoint estimateGlobalAncestry:
    log:
        OUT_DIR / "logs" / "estimateGlobalAncestry.log",
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
    conda:
        "../../envs/genomeUtils.yml"
    envmodules: *([config.get("R_module")] if config.get("R_module") else [])
    threads: 8
    resources:
        nodes=1,
        mem_mb=64000,
        runtime=240,
    input:
        labels=ancient(REF / "1000G_highcoverage" / "population.txt"),
        eigen_ref=OUT_DIR / "01-globalAncestry" / "refRefPCscores.sscore",
        eigen_sample=OUT_DIR / "01-globalAncestry" / "sampleRefPCscores.sscore",
        umap_ref=OUT_DIR / "01-globalAncestry" / "umap_ref.csv",
        umap_sample=OUT_DIR / "01-globalAncestry" / "umap_sample.csv",
        rfmix_global=OUT_DIR / "02-localAncestry" / "ancestry_full.txt" if uses_rfmix() else [],
    output:
        pos_prob=OUT_DIR / "01-globalAncestry" / "classificationProbabilities.tsv",
        sample_coords=OUT_DIR / "01-globalAncestry" / "sample_coords.tsv",
        ref_coords=OUT_DIR / "01-globalAncestry" / "ref_coords.tsv",
    params:
        dir=OUT_DIR / "01-globalAncestry",
        script=workflow.source_path("../scripts/trainPredict.R"),
        use_rfmix=uses_rfmix(),
    shell:
        """
        if [ "{params.use_rfmix}" = "True" ]; then
          rfmix_arg="--rfmix_global {input.rfmix_global}"
        else
          rfmix_arg=""
        fi
        Rscript {params.script} \
          --eigen_ref {input.eigen_ref} \
          --eigen_sample {input.eigen_sample} \
          --umap_ref {input.umap_ref} \
          --umap_sample {input.umap_sample} \
          --labels {input.labels} \
          $rfmix_arg \
          --out {params.dir} \
          --rseed $RANDOM
        """


checkpoint classifySamplesByAncestry:
    log:
        OUT_DIR / "logs" / "classifySamplesByAncestry.log",
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
    conda:
        "../../envs/genomeUtils.yml"
    envmodules: *([config.get("R_module")] if config.get("R_module") else [])
    threads: 1
    resources:
        nodes=1,
        mem_mb=16000,
        runtime=60,
    output:
        classifications=OUT_DIR / "01-globalAncestry" / "ancestry_classifications.tsv",
        keep_AFR=OUT_DIR / "01-globalAncestry" / "keep_AFR.txt",
        keep_AMR=OUT_DIR / "01-globalAncestry" / "keep_AMR.txt",
        keep_EAS=OUT_DIR / "01-globalAncestry" / "keep_EAS.txt",
        keep_EUR=OUT_DIR / "01-globalAncestry" / "keep_EUR.txt",
        keep_SAS=OUT_DIR / "01-globalAncestry" / "keep_SAS.txt",
        keep_Other=OUT_DIR / "01-globalAncestry" / "keep_Other.txt",
        ridge_plot=report(
            OUT_DIR / "01-globalAncestry" / f"classificationProbability_stacked_{ANCESTRY_MODEL}.svg",
            caption="../../report/ancestry_ridgelines.rst",
            category="Global ancestry",
        ),
        class_plot=report(
            OUT_DIR / "01-globalAncestry" / "ancestry_classification_space.svg",
            caption="../../report/ancestry_classification.rst",
            category="Global ancestry",
        ),
        density_plot=OUT_DIR / "01-globalAncestry" / f"classificationProbability_density_{ANCESTRY_MODEL}.svg",
    input:
        pos_prob=get_classification_probs,
        sample_coords=OUT_DIR / "01-globalAncestry" / "sample_coords.tsv",
        ref_coords=OUT_DIR / "01-globalAncestry" / "ref_coords.tsv",
    params:
        dir=OUT_DIR / "01-globalAncestry",
        threshold=config.get("ancestry", {}).get("threshold", 0.8),
        model=ANCESTRY_MODEL,
        script=workflow.source_path("../scripts/classify.R"),
        plot_posterior=workflow.source_path("../scripts/plotPosterior.R"),
        plot_classification=workflow.source_path("../scripts/plotClassification.R"),
        plot_density=workflow.source_path("../scripts/plotProbabilityDensity.R"),
    shell:
        """
        Rscript {params.script} \
          --out {params.dir} \
          --threshold {params.threshold} \
          --model {params.model}

        Rscript {params.plot_posterior} \
          --prob_file {input.pos_prob} \
          --out_dir {params.dir}

        Rscript {params.plot_classification} \
          --out_dir {params.dir} \
          --threshold {params.threshold} \
          --model {params.model} \
          --rf_model {params.dir}/RFpc.Rds

        Rscript {params.plot_density} \
          --prob_file {input.pos_prob} \
          --out_dir {params.dir} \
          --model {params.model}
        """


