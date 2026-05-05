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

    rule prepareSnpHeritInputs:
        container:
            "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
        conda:
            "../../envs/ancNreport.yml"
        threads: 1
        resources:
            nodes=1,
            mem_mb=8000,
            runtime=60,
        input:
            eigenvec=OUT_DIR / "01-globalAncestry" / "ref.eigenvec",
            pheno=config.get("snpHerit", {}).get("pheno"),
            covar=config.get("snpHerit", {}).get("covar"),
        output:
            config_json=OUT_DIR / "{subset}" / "03-snpHeritability" / "mash_config.json",
        params:
            npc=config.get("snpHerit", {}).get("npc", 10),
            mpheno=config.get("snpHerit", {}).get("mpheno", 1),
            method=config.get("snpHerit", {}).get("method", "AdjHE"),
            prefix=config.get("snpHerit", {}).get("prefix"),
            qcovar=config.get("snpHerit", {}).get("qcovar", None),
            covar_discrete=config.get("snpHerit", {}).get("covar_discrete", None),
            RV=config.get("snpHerit", {}).get("RV", None),
            loop_covs=config.get("snpHerit", {}).get("loop_covs", False),
            std=config.get("snpHerit", {}).get("std", False),
            k=config.get("snpHerit", {}).get("k", None),
            pheno_filter=config.get("snpHerit", {}).get("pheno_filter", None),
            covar_filter=config.get("snpHerit", {}).get("covar_filter", None),
            out_dir=lambda wildcards: OUT_DIR / wildcards.subset / "03-snpHeritability",
        run:
            import json
            import os

            os.makedirs(params.out_dir, exist_ok=True)

            # Add header to pheno if missing
            pheno_path = input.pheno
            if os.path.exists(pheno_path):
                with open(pheno_path, 'r') as f:
                    first_line = f.readline().strip()
                if not first_line.startswith(('FID', 'IID')):
                    with open(pheno_path, 'r') as f:
                        content = f.read()
                    with open(pheno_path, 'w') as f:
                        f.write('FID\tIID\tPHENO\n')
                        f.write(content)

            prefix = params.prefix
            if not prefix or str(prefix) == "None":
                prefix = os.path.join(params.out_dir, f"{wildcards.subset}_grm")
            else:
                # Replace {subset} placeholder with actual subset value
                prefix = prefix.replace("{subset}", wildcards.subset)

            config = {
                "PC": input.eigenvec,
                "prefix": prefix,
                "pheno": input.pheno,
                "covar": input.covar,
                "ids": None,
                "out": os.path.join(params.out_dir, "mash_output"),
                "npc": int(params.npc),
                "mpheno": int(params.mpheno),
                "Method": params.method
            }

            if params.qcovar is not None:
                config["qcovar"] = params.qcovar
            if params.covar_discrete is not None:
                config["covar_discrete"] = params.covar_discrete
            if params.RV and str(params.RV) != "None":
                config["RV"] = params.RV
            if params.loop_covs:
                config["loop_covars"] = True
            if params.std:
                config["std"] = True
            if params.k is not None:
                config["k"] = params.k
            if params.pheno_filter and str(params.pheno_filter) != "None":
                config["pheno_filter"] = params.pheno_filter
            if params.covar_filter and str(params.covar_filter) != "None":
                config["covar_filter"] = params.covar_filter

            with open(output.config_json, "w") as f:
                json.dump(config, f, indent=2)

    rule estimateSnpHeritability:
        container:
            "oras://ghcr.io/coffm049/gdcgenomicsqc/mash:latest"
        threads: 8
        resources:
            nodes=1,
            mem_mb=32000,
            runtime=2880,
        input:
            config_json=OUT_DIR / "{subset}" / "03-snpHeritability" / "mash_config.json",
        output:
            estimates=OUT_DIR / "{subset}" / "03-snpHeritability" / "mash_output.csv",
        log:
            OUT_DIR / "{subset}" / "03-snpHeritability" / "mash.log",
        params:
            out_dir=lambda wildcards: OUT_DIR / wildcards.subset / "03-snpHeritability",
        shell:
            """
            cd {params.out_dir}
            MASH --argfile {input.config_json} > {log} 2>&1
            """
