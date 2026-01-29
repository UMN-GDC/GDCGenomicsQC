rule plotInitial_QC:
    conda: "../../envs/ancNreport.yml"
    input:
        imiss = os.path.join(config['OUT_DIR'], "01-Initialfilter/initial.imiss"),
        lmiss = os.path.join(config['OUT_DIR'], "01-Initialfilter/initial.lmiss")
    output:
        imiss = report(os.path.join(config['OUT_DIR'], "figures/imiss.png"), caption = "../report/imiss.rst", category = "Initial QC"),
        lmiss = report(os.path.join(config['OUT_DIR'], "figures/lmiss.png"), caption = "../report/lmiss.rst", category = "Initial QC"),
    shell: """
    bash scripts/plotMissingness.R {input.imiss} {input.lmiss} {output.imiss} {output.lmiss} 
    """

