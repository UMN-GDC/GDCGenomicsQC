rule estimateAncestry:
    conda: "../../envs/ancNreport.yml"
    resources:
        nodes = 1,
        mem_mb = 64000,
        runtime = 2880,
    input:
        fam = OUT_DIR / "02-relatedness" / "standardFiltered.LDpruned.fam",
        ancestry = REF / "1000G_GRCh38" / "1000G.popu",
        eigen = OUT_DIR / "04-globalAncestry" / "merged_dataset_pca.eigenvec",
        umap = OUT_DIR / "04-globalAncestry" / "umap.csv",
    output:
        report(OUT_DIR / "04-globalAncestry" / "latentDistantRelatedness.png", caption = "../../report/PCs.rst", category = "Global ancestry"),
        OUT_DIR / "04-globalAncestry" / "latentDistantRelatedness.tsv",
    params:
        dir = OUT_DIR / "04-globalAncestry",

    shell: """
    echo "Running ancestry estimation:"

    Rscript scripts/classification.R  \
      --pc {input.eigen} \
      --umap {input.umap} \
      --labels {input.ancestry} \
      --out {params.dir} \
      --rseed $RANDOM

    """

