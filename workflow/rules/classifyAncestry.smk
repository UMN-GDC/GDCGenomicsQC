
checkpoint estimateGlobalAncestry:
    log:
        OUT_DIR / "logs" / "estimateGlobalAncestry.log",
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
    conda:
        "../../envs/genomeUtils.yml"
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
        rfmix_global=lambda wildcards: OUT_DIR / "02-localAncestry" / "ancestry_full.txt" if uses_rfmix() else [],
    output:
        pos_prob=OUT_DIR / "01-globalAncestry" / "posterior_probabilities.tsv",
        sample_coords=OUT_DIR / "01-globalAncestry" / "sample_coords.tsv",
        ref_coords=OUT_DIR / "01-globalAncestry" / "ref_coords.tsv",
        ridge_plot=report(
            OUT_DIR
            / "01-globalAncestry"
            / f"posterior_probability_stacked_{ANCESTRY_MODEL}.svg",
            caption="../../report/ancestry_ridgelines.rst",
            category="Global ancestry",
        ),
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
    threads: 1
    resources:
        nodes=1,
        mem_mb=16000,
        runtime=60,
    input:
        pos_prob=get_posterior_probs,
        sample_coords=OUT_DIR / "01-globalAncestry" / "sample_coords.tsv",
        ref_coords=OUT_DIR / "01-globalAncestry" / "ref_coords.tsv",
    output:
        classifications=OUT_DIR / "01-globalAncestry" / "ancestry_classifications.tsv",
        class_plot=report(
            OUT_DIR / "01-globalAncestry" / "ancestry_classification_space.svg",
            caption="../../report/ancestry_classification.rst",
            category="Global ancestry",
        ),
        keep_AFR=OUT_DIR / "01-globalAncestry" / "keep_AFR.txt",
        keep_AMR=OUT_DIR / "01-globalAncestry" / "keep_AMR.txt",
        keep_EAS=OUT_DIR / "01-globalAncestry" / "keep_EAS.txt",
        keep_EUR=OUT_DIR / "01-globalAncestry" / "keep_EUR.txt",
        keep_SAS=OUT_DIR / "01-globalAncestry" / "keep_SAS.txt",
        keep_Other=OUT_DIR / "01-globalAncestry" / "keep_Other.txt",
    params:
        dir=OUT_DIR / "01-globalAncestry",
        threshold=config.get("ancestry", {}).get("threshold", 0.8),
        model=ANCESTRY_MODEL,
        script=workflow.source_path("../scripts/classify.R"),
    shell:
        """
        Rscript {params.script} \
          --out {params.dir} \
          --threshold {params.threshold} \
          --model {params.model}
        """
