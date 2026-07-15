SNP_HERIT_CONFIG = config.get("snpHerit", {})
SNP_HERIT_ACTIVE = bool(SNP_HERIT_CONFIG.get("pheno"))
SNP_HERIT_OUT = SNP_HERIT_CONFIG.get("out")

if SNP_HERIT_CONFIG:
    if SNP_HERIT_CONFIG.get("covar") and not SNP_HERIT_CONFIG.get("pheno"):
        raise ValueError("snpHerit.pheno must be specified in config when covar is specified")
    valid_methods = ["AdjHE", "AdjHE_fixed", "AdjHE_mixed", "AdjHE_random", "GCTA", "PredLMM", "SWD", "Combat", "Covbat"]
    method = SNP_HERIT_CONFIG.get("method", "AdjHE")
    if method not in valid_methods:
        raise ValueError(f"snpHerit.method must be one of {valid_methods}")

import json


def _snp_herit_out_dir(w):
    if SNP_HERIT_OUT:
        return Path(SNP_HERIT_OUT)
    return OUT_DIR / w.subset / "03-snpHeritability"


def _mash_config(prefix, pheno, out, npc, mpheno, eigenvec,
                 covar=None, covar_discrete=None, qcovar=None,
                 pheno_filter=None, covar_filter=None,
                 loop_covars=False, random_groups=None, Naive=False,
                 std=None, k=None, RV=None):
    cfg = {
        "prefix": str(prefix),
        "pheno": [str(p) for p in (pheno if isinstance(pheno, list) else [pheno])],
        "out": str(out),
        "npc": [int(n) for n in (npc if isinstance(npc, list) else [npc])],
        "mpheno": [str(m) for m in (mpheno if isinstance(mpheno, list) else [mpheno])],
        "Method": SNP_HERIT_CONFIG.get("method", "AdjHE"),
        "iid_col": SNP_HERIT_CONFIG.get("iid_col", "IID"),
        "fid_col": SNP_HERIT_CONFIG.get("fid_col", "FID"),
        "PC": str(eigenvec),
        "loop_covars": loop_covars,
        "random_groups": RV if RV else None,
        "Naive": Naive,
    }
    if covar:
        cfg["covar"] = [str(c) for c in covar] if isinstance(covar, list) else str(covar)
    if covar_discrete:
        cfg["covar_discrete"] = covar_discrete
    if qcovar:
        cfg["qcovar"] = qcovar
    if pheno_filter:
        cfg["pheno_filter"] = pheno_filter
    if covar_filter:
        cfg["covar_filter"] = covar_filter
    if std:
        cfg["std"] = std
    if k is not None:
        cfg["k"] = k
    if RV:
        cfg["RV"] = RV
    return json.dumps(cfg, indent=2)


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
            grm_bin=OUT_DIR / "{subset}" / "f1.b38.ldpruned.unrelated.grm.bin",
            grm_id=OUT_DIR / "{subset}" / "f1.b38.ldpruned.unrelated.grm.id",
            grm_Nbin=OUT_DIR / "{subset}" / "f1.b38.ldpruned.unrelated.grm.N.bin",
            eigenvec=OUT_DIR / "{subset}" / "internal_pca_plink2.eigenvec",
        output:
            estimates=lambda w: _snp_herit_out_dir(w) / "mash_output.csv",
        params:
            argfile=lambda w: _snp_herit_out_dir(w) / "mash_output.json",
            mash_config=lambda w: _mash_config(
                prefix=OUT_DIR / w.subset / "f1.b38.ldpruned.unrelated",
                pheno=SNP_HERIT_CONFIG["pheno"],
                out=_snp_herit_out_dir(w) / "mash_output",
                npc=SNP_HERIT_CONFIG.get("npc", 10),
                mpheno=SNP_HERIT_CONFIG.get("mpheno", 1),
                eigenvec=OUT_DIR / w.subset / "internal_pca_plink2.eigenvec",
                covar=SNP_HERIT_CONFIG.get("covar"),
                qcovar=SNP_HERIT_CONFIG.get("qcovar"),
                covar_discrete=SNP_HERIT_CONFIG.get("covar_discrete"),
                pheno_filter=SNP_HERIT_CONFIG.get("pheno_filter"),
                covar_filter=SNP_HERIT_CONFIG.get("covar_filter"),
                loop_covars=SNP_HERIT_CONFIG.get("loop_covars", False),
                random_groups=SNP_HERIT_CONFIG.get("RV"),
                Naive=SNP_HERIT_CONFIG.get("Naive", False),
                std=SNP_HERIT_CONFIG.get("std"),
                k=SNP_HERIT_CONFIG.get("k"),
                RV=SNP_HERIT_CONFIG.get("RV"),
            ),
        shell:
            """
            mkdir -p "$(dirname {output.estimates})"
            cat > {params.argfile} << 'EOF'
{params.mash_config}
EOF
            MASH --argfile {params.argfile}
            """

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
            grm_bin=OUT_DIR / "{subset}" / "simulations" / "{sim_name}" / "simulated.grm.bin",
            grm_id=OUT_DIR / "{subset}" / "simulations" / "{sim_name}" / "simulated.grm.id",
            grm_Nbin=OUT_DIR / "{subset}" / "simulations" / "{sim_name}" / "simulated.grm.N.bin",
            eigenvec=OUT_DIR / "{subset}" / "simulations" / "{sim_name}" / "simulated.eigenvec",
            pheno=OUT_DIR / "{subset}" / "simulations" / "{sim_name}" / "simulated_pheno1.pheno",
        output:
            estimates=OUT_DIR / "{subset}" / "simulations" / "{sim_name}" / "herit.csv",
        params:
            argfile=lambda w: OUT_DIR / w.subset / "simulations" / w.sim_name / "mash_config.json",
            mash_config=lambda w: _mash_config(
                prefix=OUT_DIR / w.subset / "simulations" / w.sim_name / "simulated",
                pheno=[OUT_DIR / w.subset / "simulations" / w.sim_name / "simulated_pheno1.pheno"],
                out=OUT_DIR / w.subset / "simulations" / w.sim_name / "herit",
                npc=SNP_HERIT_CONFIG.get("npc", 10),
                mpheno=SNP_HERIT_CONFIG.get("mpheno", 1),
                eigenvec=OUT_DIR / w.subset / "simulations" / w.sim_name / "simulated.eigenvec",
                covar=SNP_HERIT_CONFIG.get("covar"),
                qcovar=SNP_HERIT_CONFIG.get("qcovar"),
                covar_discrete=SNP_HERIT_CONFIG.get("covar_discrete"),
                pheno_filter=SNP_HERIT_CONFIG.get("pheno_filter"),
                covar_filter=SNP_HERIT_CONFIG.get("covar_filter"),
                loop_covars=SNP_HERIT_CONFIG.get("loop_covars", False),
                random_groups=SNP_HERIT_CONFIG.get("RV"),
                Naive=SNP_HERIT_CONFIG.get("Naive", False),
                std=SNP_HERIT_CONFIG.get("std"),
                k=SNP_HERIT_CONFIG.get("k"),
                RV=SNP_HERIT_CONFIG.get("RV"),
            ),
        shell:
            """
            mkdir -p "$(dirname {output.estimates})"
            cat > {params.argfile} << 'EOF'
{params.mash_config}
EOF
            MASH --argfile {params.argfile}
            """
