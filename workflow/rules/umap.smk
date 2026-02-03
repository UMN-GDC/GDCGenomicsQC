rule UMAP:
    conda: "../../envs/ancNreport.yml"
    threads: 8
    resources:
        nodes = 1,
        mem_mb = 64000,
        runtime = 2880,
        slurm_extra = "'--job-name=UMAP'"
    input:
        eigen = f"{config['OUT_DIR']}/04-globalAncestry/merged_dataset_pca.eigenvec",
    output:
        os.path.join(config['OUT_DIR'], "04-globalAncestry/umap.csv"),
    params:
        npc = 10,
        neighbors = 50,
        ncoords = 2,
    shell: """
    echo "Running UMAP:"

    Rscript scripts/Umap.R --eigens {input.eigen} --out {output} \
      --npc {params.npc} --neighbors {params.neighbors} \
      --threads {threads} \
      --ncoords {params.ncoords} \
      --seed $RANDOM

    """

