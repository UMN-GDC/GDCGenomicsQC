rule plotPCs:
    conda: "../../envs/ancNreport.yml"
    input:
        eigen = f"{config['OUT_DIR']}/04-globalAncestry/{subset}.eigenvec",
        fam = f"{config['OUT_DIR']}/04-globalAncestry/merged_common_bi.fam",
        ancestry = f"{config['OUT_DIR']}/04-globalAncestry/data.txt",
    output:
        report(os.path.join(config['OUT_DIR'], "figures/{subset}PCs.png"), caption = "../report/lmiss.rst", category = "Initial QC"),
        popu = lambda wildcards: get_input_by_stage(wildcards) + ".popu",
    shell: """


    Rscript scripts/plot_pca.R {input.eigen} {input.fam} {input.popu} {output}

    """

