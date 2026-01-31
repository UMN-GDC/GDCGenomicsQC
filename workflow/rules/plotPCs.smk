rule plotInitial_QC:
    conda: "../../envs/ancNreport.yml"
    input:
        eigen = f"{config['OUT_DIR']}/04-globalAncestry/merged_dataset_pca.eigenvec",
    output:
        report(os.path.join(config['OUT_DIR'], "figures/PCs.png"), caption = "../report/lmiss.rst", category = "Initial QC"),
    shell: """
    Rscript scripts/plotMissingness.R {input.imiss} {input.lmiss} {output.imiss} {output.lmiss} 
    """

