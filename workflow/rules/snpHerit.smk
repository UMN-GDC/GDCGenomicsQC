SNP_HERIT_CONFIG = config.get("snpHerit", {})
SNP_HERIT_ACTIVE = bool(SNP_HERIT_CONFIG.get("pheno"))

if SNP_HERIT_CONFIG:
    if SNP_HERIT_CONFIG.get("covar") and not SNP_HERIT_CONFIG.get("pheno"):
        raise ValueError("snpHerit.pheno must be specified in config when covar is specified")
    valid_methods = ["AdjHE", "GCTA", "PredLMM", "SWD", "Combat", "Covbat"]
    method = SNP_HERIT_CONFIG.get("method", "AdjHE")
    if method not in valid_methods:
        raise ValueError(f"snpHerit.method must be one of {valid_methods}")

if SNP_HERIT_ACTIVE:

    def get_snpHerit_pca_input(wildcards):
        pca_input = SNP_HERIT_CONFIG.get("pca_input")
        if pca_input:
            return {"pca": pca_input}
        return {}

    def get_snpHerit_covar_input(wildcards):
        covar = SNP_HERIT_CONFIG.get("covar")
        if covar:
            return {"covar": covar}
        return {}

    rule prepareSnpHeritArgfile:
        input:
            lambda wildcards: get_snpHerit_pca_input(wildcards),
            lambda wildcards: get_snpHerit_covar_input(wildcards),
            pheno=config.get("snpHerit", {}).get("pheno"),
        output:
            argfile=OUT_DIR / "{subset}" / "03-snpHeritability" / "mash_argfile.json",
        params:
            out_dir=lambda wildcards: OUT_DIR / wildcards.subset / "03-snpHeritability",
            mash_out=lambda wildcards: (
                str(Path(SNP_HERIT_CONFIG.get("out")).with_suffix(""))
                if SNP_HERIT_CONFIG.get("out")
                else str(OUT_DIR / wildcards.subset / "03-snpHeritability") + "/mash_output"
            ),
            prefix=lambda wildcards: (
                SNP_HERIT_CONFIG.get("grm_prefix")
                if SNP_HERIT_CONFIG.get("grm_prefix")
                else OUT_DIR / wildcards.subset / f"{wildcards.subset}_grm"
            ),
            pca_input=SNP_HERIT_CONFIG.get("pca_input", None),
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

            mash_config = {
                "prefix": str(params.prefix),
                "pheno": str(input.pheno),
                "out": params.mash_out,
                "npc": int(params.npc),
                "mpheno": int(params.mpheno),
                "Method": params.method,
                "iid_col": params.iid_col,
                "fid_col": params.fid_col,
            }

            if "covar" in input:
                mash_config["covar"] = str(input.covar)

            if params.pca_input:
                pca_path = str(params.pca_input)
                if pca_path.endswith(".RDS"):
                    mash_config["PC"] = pca_path
                elif pca_path.endswith(".eigenvec") or pca_path.endswith(".eigenvec.txt"):
                    mash_config["eigenvec"] = pca_path

            if params.qcovar:
                mash_config["qcovar"] = params.qcovar
            if params.covar_discrete:
                mash_config["covar_discrete"] = params.covar_discrete
            if params.pheno_filter:
                mash_config["pheno_filter"] = params.pheno_filter
            if params.covar_filter:
                mash_config["covar_filter"] = params.covar_filter

            with open(output.argfile, "w") as f:
                json.dump(mash_config, f, indent=2)

    rule estimateSnpHeritability:
        log:
            OUT_DIR / "{subset}" / "03-snpHeritability" / "mash.log",
        container:
            "oras://ghcr.io/coffm049/gdcgenomicsqc/mash:v1"
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
            mkdir -p {output.estimates.parent}
            MASH estimate --argfile {input.argfile} > {log} 2>&1
            """
