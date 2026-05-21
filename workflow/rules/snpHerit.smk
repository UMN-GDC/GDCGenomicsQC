SNP_HERIT_CONFIG = config.get("snpHerit", {})
SNP_HERIT_ACTIVE = bool(SNP_HERIT_CONFIG.get("pheno"))

if SNP_HERIT_CONFIG:
    if SNP_HERIT_CONFIG.get("covar") and not SNP_HERIT_CONFIG.get("pheno"):
        raise ValueError("snpHerit.pheno must be specified in config when covar is specified")
    valid_methods = ["AdjHE", "AdjHE_fixed", "AdjHE_mixed", "AdjHE_random", "GCTA", "PredLMM", "SWD", "Combat", "Covbat"]
    method = SNP_HERIT_CONFIG.get("method", "AdjHE")
    if method not in valid_methods:
        raise ValueError(f"snpHerit.method must be one of {valid_methods}")

if SNP_HERIT_ACTIVE:

    rule prepareSnpHeritArgfile:
        input:
            pca=lambda wildcards: [SNP_HERIT_CONFIG["pca_input"]] if SNP_HERIT_CONFIG.get("pca_input") else [],
            covar=lambda wildcards: list(SNP_HERIT_CONFIG["covar"]) if SNP_HERIT_CONFIG.get("covar") else [],
            pheno=lambda wildcards: list(SNP_HERIT_CONFIG["pheno"]) if SNP_HERIT_CONFIG.get("pheno") else [],
        output:
            argfile=Path(str(SNP_HERIT_CONFIG["out"]).replace(".csv", ".json")).resolve(),
        run:
            import json
            import os

            pheno_val = input.pheno
            if isinstance(pheno_val, list):
                pheno_list = [str(p) for p in pheno_val]
            else:
                pheno_list = [str(pheno_val)]

            mpheno_val = SNP_HERIT_CONFIG.get("mpheno", 1)
            if isinstance(mpheno_val, list):
                mpheno_list = [str(m) for m in mpheno_val]
            else:
                mpheno_list = [str(mpheno_val)]

            mash_config = {
                "prefix": str(SNP_HERIT_CONFIG.get("grm_prefix")),
                "pheno": pheno_list,
                "out": str(Path(str(SNP_HERIT_CONFIG["out"])).with_suffix("")),
                "npc": [int(SNP_HERIT_CONFIG.get("npc", 10))] if not isinstance(SNP_HERIT_CONFIG.get("npc", 10), list) else [int(n) for n in SNP_HERIT_CONFIG["npc"]],
                "mpheno": mpheno_list,
                "Method": SNP_HERIT_CONFIG.get("method", "AdjHE"),
                "iid_col": SNP_HERIT_CONFIG.get("iid_col", "IID"),
                "fid_col": SNP_HERIT_CONFIG.get("fid_col", "FID"),
            }

            covar_from_config = SNP_HERIT_CONFIG.get("covar")
            if covar_from_config:
                if isinstance(covar_from_config, list):
                    mash_config["covar"] = [str(f) for f in covar_from_config]
                else:
                    mash_config["covar"] = str(covar_from_config)

            if SNP_HERIT_CONFIG.get("pca_input"):
                mash_config["PC"] = str(SNP_HERIT_CONFIG["pca_input"])

            mash_config["qcovar"] = SNP_HERIT_CONFIG.get("qcovar")
            mash_config["covar_discrete"] = SNP_HERIT_CONFIG.get("covar_discrete")
            if SNP_HERIT_CONFIG.get("pheno_filter"):
                mash_config["pheno_filter"] = SNP_HERIT_CONFIG["pheno_filter"]
            if SNP_HERIT_CONFIG.get("covar_filter"):
                mash_config["covar_filter"] = SNP_HERIT_CONFIG["covar_filter"]

            mash_config["loop_covars"] = SNP_HERIT_CONFIG.get("loop_covars", False)
            mash_config["random_groups"] = SNP_HERIT_CONFIG.get("RV", None) if SNP_HERIT_CONFIG.get("RV") else None
            mash_config["Naive"] = SNP_HERIT_CONFIG.get("Naive", False)

            if SNP_HERIT_CONFIG.get("std"):
                mash_config["std"] = SNP_HERIT_CONFIG["std"]
            if SNP_HERIT_CONFIG.get("k") is not None:
                mash_config["k"] = SNP_HERIT_CONFIG["k"]
            if SNP_HERIT_CONFIG.get("RV"):
                mash_config["RV"] = SNP_HERIT_CONFIG["RV"]

            os.makedirs(os.path.dirname(str(output.argfile)), exist_ok=True)
            with open(str(output.argfile), "w") as f:
                json.dump(mash_config, f, indent=2)

    rule estimateSnpHeritability:
        log:
            str(Path(str(SNP_HERIT_CONFIG["out"]).replace(".csv", ".log")).resolve()),
        conda:
            "../../envs/mash.yml"
        container:
            "oras://ghcr.io/coffm049/gdcgenomicsqc/mash:v1"
        envmodules:
            mod("R")
        threads: 8
        resources:
            nodes=1,
            mem_mb=32000,
            runtime=2880,
        input:
            argfile=rules.prepareSnpHeritArgfile.output.argfile,
        output:
            estimates=str(Path(str(SNP_HERIT_CONFIG["out"])).resolve()),
        shell:
            """
            MASH --argfile {input.argfile} > {log} 2>&1
            """
