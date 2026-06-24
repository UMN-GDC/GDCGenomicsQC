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

    rule estimateSnpHeritability:
        conda:
            "../../envs/mash.yml"
        container:
            "oras://ghcr.io/coffm049/gdcgenomicsqc/mash:v1"
        envmodules: *([config.get("R_module")] if config.get("R_module") else [])
        threads: 8
        resources:
            nodes=1,
            mem_mb=32000,
            runtime=720,
        input:
            grm_bin=OUT_DIR / "{ancestry}" / "f1.b38.ldpruned.unrelated.grm.bin",
            grm_id=OUT_DIR / "{ancestry}" / "f1.b38.ldpruned.unrelated.grm.id",
            grm_Nbin=OUT_DIR / "{ancestry}" / "f1.b38.ldpruned.unrelated.grm.N.bin",
            eigenvec=OUT_DIR / "{ancestry}" / "internal_pca_plink2.eigenvec",
        output:
            estimates=OUT_DIR / "{ancestry}" / "03-snpHeritability" / "mash_output.csv",
        params:
            grm_prefix=OUT_DIR / "{ancestry}" / "f1.b38.ldpruned.unrelated",
            out_prefix=OUT_DIR / "{ancestry}" / "03-snpHeritability" / "mash_output",
        run:
            import json
            import os

            pheno_val = SNP_HERIT_CONFIG["pheno"]
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
                "prefix": str(params.grm_prefix),
                "pheno": pheno_list,
                "out": str(params.out_prefix),
                "npc": [int(SNP_HERIT_CONFIG.get("npc", 10))] if not isinstance(SNP_HERIT_CONFIG.get("npc", 10), list) else [int(n) for n in SNP_HERIT_CONFIG["npc"]],
                "mpheno": mpheno_list,
                "Method": SNP_HERIT_CONFIG.get("method", "AdjHE"),
                "iid_col": SNP_HERIT_CONFIG.get("iid_col", "IID"),
                "fid_col": SNP_HERIT_CONFIG.get("fid_col", "FID"),
                "PC": str(input.eigenvec),
            }

            covar_from_config = SNP_HERIT_CONFIG.get("covar")
            if covar_from_config:
                if isinstance(covar_from_config, list):
                    mash_config["covar"] = [str(f) for f in covar_from_config]
                else:
                    mash_config["covar"] = str(covar_from_config)

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

            os.makedirs(os.path.dirname(str(output.estimates)), exist_ok=True)
            argfile = str(output.estimates).replace(".csv", ".json")
            with open(argfile, "w") as f:
                json.dump(mash_config, f, indent=2)

            import subprocess
            subprocess.run(["MASH", "--argfile", argfile], check=True)


SIM_CFG = config.get("phenotypeSimulation", {})
if SIM_CFG.get("enabled", False):
    if not SIM_CFG.get("ancestries"):
        raise ValueError("phenotypeSimulation.ancestries must be specified when enabled")

    def get_sim_cfg(sim_name):
        for s in SIM_CFG.get("simulations", []):
            if s.get("name") == sim_name:
                return s
        return {}

    rule estimateSnpHeritabilitySimulated:
        conda:
            "../../envs/mash.yml"
        container:
            "oras://ghcr.io/coffm049/gdcgenomicsqc/mash:v1"
        envmodules: *([config.get("R_module")] if config.get("R_module") else [])
        threads: 8
        resources:
            nodes=1,
            mem_mb=32000,
            runtime=720,
        input:
            grm_bin=OUT_DIR / "{ancestry}" / "simulations" / "{sim_name}" / "simulated.grm.bin",
            grm_id=OUT_DIR / "{ancestry}" / "simulations" / "{sim_name}" / "simulated.grm.id",
            grm_Nbin=OUT_DIR / "{ancestry}" / "simulations" / "{sim_name}" / "simulated.grm.N.bin",
            eigenvec=OUT_DIR / "{ancestry}" / "simulations" / "{sim_name}" / "simulated.eigenvec",
            pheno=OUT_DIR / "{ancestry}" / "simulations" / "{sim_name}" / "simulated_pheno1.pheno",
        output:
            estimates=OUT_DIR / "{ancestry}" / "simulations" / "{sim_name}" / "herit.csv",
        params:
            grm_prefix=OUT_DIR / "{ancestry}" / "simulations" / "{sim_name}" / "simulated",
            out_prefix=OUT_DIR / "{ancestry}" / "simulations" / "{sim_name}" / "herit",
        run:
            import json
            import os

            mpheno_val = SNP_HERIT_CONFIG.get("mpheno", 1)
            if isinstance(mpheno_val, list):
                mpheno_list = [str(m) for m in mpheno_val]
            else:
                mpheno_list = [str(mpheno_val)]

            mash_config = {
                "prefix": str(params.grm_prefix),
                "pheno": [str(input.pheno)],
                "out": str(params.out_prefix),
                "npc": [int(SNP_HERIT_CONFIG.get("npc", 10))] if not isinstance(SNP_HERIT_CONFIG.get("npc", 10), list) else [int(n) for n in SNP_HERIT_CONFIG["npc"]],
                "mpheno": mpheno_list,
                "Method": SNP_HERIT_CONFIG.get("method", "AdjHE"),
                "iid_col": SNP_HERIT_CONFIG.get("iid_col", "IID"),
                "fid_col": SNP_HERIT_CONFIG.get("fid_col", "FID"),
                "PC": str(input.eigenvec),
            }

            covar_from_config = SNP_HERIT_CONFIG.get("covar")
            if covar_from_config:
                if isinstance(covar_from_config, list):
                    mash_config["covar"] = [str(f) for f in covar_from_config]
                else:
                    mash_config["covar"] = str(covar_from_config)

            mash_config["qcovar"] = SNP_HERIT_CONFIG.get("qcovar")
            mash_config["covar_discrete"] = SNP_HERIT_CONFIG.get("covar_discrete")

            mash_config["loop_covars"] = SNP_HERIT_CONFIG.get("loop_covars", False)
            mash_config["random_groups"] = SNP_HERIT_CONFIG.get("RV", None) if SNP_HERIT_CONFIG.get("RV") else None
            mash_config["Naive"] = SNP_HERIT_CONFIG.get("Naive", False)

            if SNP_HERIT_CONFIG.get("std"):
                mash_config["std"] = SNP_HERIT_CONFIG["std"]
            if SNP_HERIT_CONFIG.get("k") is not None:
                mash_config["k"] = SNP_HERIT_CONFIG["k"]
            if SNP_HERIT_CONFIG.get("RV"):
                mash_config["RV"] = SNP_HERIT_CONFIG["RV"]

            os.makedirs(os.path.dirname(str(output.estimates)), exist_ok=True)
            argfile = str(output.estimates).replace(".csv", ".json")
            with open(argfile, "w") as f:
                json.dump(mash_config, f, indent=2)

            import subprocess
            subprocess.run(["MASH", "--argfile", argfile], check=True)
