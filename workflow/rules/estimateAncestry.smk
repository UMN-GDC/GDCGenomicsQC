rule estimateAncestry:
    container: "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
    conda: "../../envs/ancNreport.yml"
    resources:
        nodes = 1,
        mem_mb = 64000,
        runtime = 2880,
    input:
        labels = REF / "1000G_GRCh38"/ "1000G.GRCh38.popu",
        eigen_ref = OUT_DIR / "01-globalAncestry" / "refRefPCscores.sscore",
        eigen_sample = OUT_DIR / "01-globalAncestry" / "sampleRefPCscores.sscore",
        umap_ref = OUT_DIR / "01-globalAncestry" / "umap_ref.csv",
        umap_sample = OUT_DIR / "01-globalAncestry" / "umap_sample.csv",
    output:
        report(OUT_DIR / "01-globalAncestry" / "PC_referenceSpace.svg", caption = "../../report/PCs.rst", category = "Global ancestry"),
        report(OUT_DIR / "01-globalAncestry" / "UMAP_referenceSpace.svg", caption = "../../report/PCs.rst", category = "Global ancestry"),
        OUT_DIR / "01-globalAncestry" / "latentDistantRelatedness.tsv",
    params:
        dir = OUT_DIR / "01-globalAncestry",


    shell: """
    echo "Running ancestry estimation:"

    Rscript scripts/classification.R  \
      --eigen_ref {input.eigen_ref} \
      --eigen_sample {input.eigen_sample} \
      --umap_ref {input.umap_ref} \
      --umap_sample {input.umap_sample} \
      --labels {input.labels} \
      --out {params.dir} \
      --rseed $RANDOM

    """

