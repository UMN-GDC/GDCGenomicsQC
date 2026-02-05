rule estimateAncestry:
    conda: "../../envs/ancNreport.yml"
    threads: 8
    resources:
        nodes = 1,
        mem_mb = 64000,
        runtime = 2880,
        slurm_extra = "'--job-name=estimateAncestry'"
    input:
        fam = os.path.join(config['OUT_DIR'], "02-relatedness/standardFiltered.LDpruned.fam"),
        ancestry = os.path.join(config['REF'], "rfmix_ref/super_population_map_file.txt"),
        eigen = os.path.join(config['OUT_DIR'], "04-globalAncestry/merged_dataset_pca.eigenvec"),
        umap = os.path.join(config['OUT_DIR'], "04-globalAncestry/umap.csv"),
    output:
        report(os.path.join(config['OUT_DIR'], "04-globalAncestry/latentDistantRelatedness.png"), caption = "../../report/PCs.rst", category = "Global ancestry"),
        os.path.join(config['OUT_DIR'], "04-globalAncestry/latentDistantRelatedness.csv"),
    params:
        dir = os.path.join(config['OUT_DIR'], "04-globalAncestry"),

    shell: """
    echo "Running ancestry estimation:"

    Rscript scripts/classification.R  \
      --pc {input.eigen} \
      --umap {input.umap} \
      --labels {input.ancestry} \
      --out {params.dir} \
      --rseed $RANDOM

    """

