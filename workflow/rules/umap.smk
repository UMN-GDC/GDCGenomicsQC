rule UMAP:
    log:
        OUT_DIR / "logs" / "UMAP.log",
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
    conda:
        "../../envs/ancNreport.yml"
    threads: 8
    resources:
        nodes=1,
        mem_mb=64000,
        runtime=2880,
    input:
        eigen=OUT_DIR / "01-globalAncestry" / "refRefPCscores.sscore",
        sample=OUT_DIR / "01-globalAncestry" / "sampleRefPCscores.sscore",
    output:
        OUT_DIR / "01-globalAncestry" / "umap_sample.csv",
        OUT_DIR / "01-globalAncestry" / "umap_ref.csv",
    params:
        npc=10,
        neighbors=50,
        ncoords=2,
        outputPrefix=OUT_DIR / "01-globalAncestry" / "umap",
    shell:
        """
    echo "Running UMAP:"

    Rscript scripts/Umap.R --eigens {input.eigen} --out {params.outputPrefix} \
      --npc {params.npc} --neighbors {params.neighbors} \
      --sample {input.sample} \
      --threads {threads} \
      --ncoords {params.ncoords} \
      --seed $RANDOM

    """
