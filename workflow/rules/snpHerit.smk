SNP_HERIT_CONFIG = config.get("snpHerit", {})
SNP_HERIT_ACTIVE = bool(SNP_HERIT_CONFIG.get("pheno") and SNP_HERIT_CONFIG.get("covar"))

if SNP_HERIT_CONFIG:
    if SNP_HERIT_CONFIG.get("pheno") and not SNP_HERIT_CONFIG.get("covar"):
        raise ValueError("snpHerit.covar must be specified in config when pheno is specified")
    if SNP_HERIT_CONFIG.get("covar") and not SNP_HERIT_CONFIG.get("pheno"):
        raise ValueError("snpHerit.pheno must be specified in config when covar is specified")
    valid_methods = ["AdjHE", "GCTA", "PredLMM", "SWD", "Combat", "Covbat"]
    method = SNP_HERIT_CONFIG.get("method", "AdjHE")
    if method not in valid_methods:
        raise ValueError(f"snpHerit.method must be one of {valid_methods}")

if SNP_HERIT_ACTIVE:

    rule prepareSnpHeritArgfile:
        input:
            pcaobj=OUT_DIR / "{subset}" / "pcair_pcaobj.RDS",
            pheno=config.get("snpHerit", {}).get("pheno"),
            covar=config.get("snpHerit", {}).get("covar"),
        output:
            argfile=OUT_DIR / "{subset}" / "03-snpHeritability" / "mash_argfile.json",
        params:
            out_dir=lambda wildcards: OUT_DIR / wildcards.subset / "03-snpHeritability",
            prefix=lambda wildcards: (
                SNP_HERIT_CONFIG.get("grm_prefix")
                if SNP_HERIT_CONFIG.get("grm_prefix")
                else OUT_DIR / wildcards.subset / f"{wildcards.subset}_grm"
            ),
            eigenvec=SNP_HERIT_CONFIG.get("eigenvec", None),
            npc=SNP_HERIT_CONFIG.get("npc", 10),
            mpheno=SNP_HERIT_CONFIG.get("mpheno", 1),
            method=SNP_HERIT_CONFIG.get("method", "AdjHE"),
            qcovar=SNP_HERIT_CONFIG.get("qcovar", None),
            covar_discrete=SNP_HERIT_CONFIG.get("covar_discrete", None),
            iid_col=SNP_HERIT_CONFIG.get("iid_col", "IID"),
            fid_col=SNP_HERIT_CONFIG.get("fid_col", "FID"),
            pheno_filter=SNP_HERIT_CONFIG.get("pheno_filter", None),
            covar_filter=SNP_HERIT_CONFIG.get("covar_filter", None),
        run:
            import json
            import os

            os.makedirs(params.out_dir, exist_ok=True)

            config = {
                "prefix": str(params.prefix),
                "PC": str(input.pcaobj),
                "pheno": str(input.pheno),
                "covar": str(input.covar),
                "out": str(params.out_dir) + "/mash_output",
                "npc": int(params.npc),
                "mpheno": int(params.mpheno),
                "Method": params.method,
                "iid_col": params.iid_col,
                "fid_col": params.fid_col,
            }

            if params.qcovar:
                config["qcovar"] = params.qcovar
            if params.covar_discrete:
                config["covar_discrete"] = params.covar_discrete
            if params.pheno_filter:
                config["pheno_filter"] = params.pheno_filter
            if params.covar_filter:
                config["covar_filter"] = params.covar_filter
            if params.eigenvec:
                config["eigenvec"] = params.eigenvec

            with open(output.argfile, "w") as f:
                json.dump(config, f, indent=2)

    rule estimateSnpHeritability:
        log:
            OUT_DIR / "{subset}" / "03-snpHeritability" / "mash.log",
        container:
            "oras://ghcr.io/coffm049/gdcgenomicsqc/mash:latest"
        threads: 8
        resources:
            nodes=1,
            mem_mb=32000,
            runtime=2880,
        input:
            argfile=rules.prepareSnpHeritArgfile.output.argfile,
        output:
            estimates=OUT_DIR / "{subset}" / "03-snpHeritability" / "mash_output.csv",
        shell:
            """
            mkdir -p {output.estimates".parent"}
            MASH estimate --argfile {input.argfile} > {log} 2>&1
            """
