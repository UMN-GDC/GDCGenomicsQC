rule snpHerit:
    container: "oras://ghcr.io/coffm049/gdcgenomicsqc/mash:latest"
    threads: 1
    resources:
        nodes = 1,
        mem_mb = 32000,
        runtime = 2880,
    input:
        grm = OUT_DIR / "full" / "initialFilter.grm.bin",
        grmid = OUT_DIR / "full" / "initialFilter.grm.id",
        grmN = OUT_DIR / "full" / "initialFilter.grm.N",
        eigen = OUT_DIR / "full" / "initialFilter.eigenvec",
        covar = OUT_DIR / "full" / "initialFilter.bim",
        pheno = OUT_DIR / "full" / "initialFilter.fam",
    output:
        estimates = OUT_DIR / "03-snpHeritability" / config['snpHerit']['out'],
    params:
        method = config['snpHerit']["method"],
        randomEffects = config['snpHerit']["randomEffects"],
        grm_prefix = lambda wildcards, input: input.grm[:-8],
        npc = config["snpHerit"]["npc"],
        covars = config["snpHerit"]["covars"],
        mpheno = config["snpHerit"]["mpheno"],
        mpheno = config["snpHerit"]["method"]
    shell: 
        """
        echo "Estimating SNP heritability: "
        mkdir -p {output.tempDir}

        MASH --PC {input.eigen} \
          --covar {input.covar} \
          --prefix {params.grm_prefix} \
          --pheno {input.pheno} \
          --out {output.estimates} \
          --npc {params.npc} \
          --mpheno {params.mpheno} \
          --covars {params.covars} \
          --Method {params.method}
        """

