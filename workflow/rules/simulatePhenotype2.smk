import json

SIM_CFG = config.get("phenotypeSimulation", {})
SIM_ANCESTRIES = SIM_CFG.get("ancestries", [])
SIM_INPUT_PREFIXES = SIM_CFG.get("input_prefixes") or {}

if not SIM_ANCESTRIES:
    raise ValueError("phenotypeSimulation.ancestries must be specified")

def get_sim_input_prefix(anc):
    if anc in SIM_INPUT_PREFIXES:
        return SIM_INPUT_PREFIXES[anc]
    return str(OUT_DIR / anc / "f1")

def get_sim_cfg(sim_name):
    for s in SIM_CFG.get("simulations", []):
        if s.get("name") == sim_name:
            return s
    return {}

_sim_inputs = {}
for anc in SIM_ANCESTRIES:
    p = get_sim_input_prefix(anc)
    _sim_inputs[f"pgen_{anc}"] = Path(p + ".pgen")
    _sim_inputs[f"pvar_{anc}"] = Path(p + ".pvar")
    _sim_inputs[f"psam_{anc}"] = Path(p + ".psam")

_sim_outputs = []
for anc in SIM_ANCESTRIES:
    _sim_outputs.append(OUT_DIR / anc / "simulations" / "{sim_name}" / "simulated.bed")
    _sim_outputs.append(OUT_DIR / anc / "simulations" / "{sim_name}" / "simulated.bim")
    _sim_outputs.append(OUT_DIR / anc / "simulations" / "{sim_name}" / "simulated.fam")

rule simulatePhenotypes:
    log:
        OUT_DIR / "logs" / "simulatePhenotypes_{sim_name}.log",
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/phenotypesim:v1"
    conda:
        "../../envs/phenotypeSim.yml"
    threads: 8
    resources:
        nodes=1,
        mem_mb=32000,
        runtime=240,
    input:
        **_sim_inputs,
    output:
        _sim_outputs,
    run:
        anc_names = SIM_ANCESTRIES
        pgen_prefixes = [get_sim_input_prefix(anc) for anc in anc_names]
        out_dirs = [str(OUT_DIR / anc / "simulations" / wildcards.sim_name) for anc in anc_names]
        sim = get_sim_cfg(wildcards.sim_name)

        import subprocess
        subprocess.run([
            "Rscript", str(SCRIPTS_DIR / "runPhenotypeSimulation.R"),
            "--anc-names", ",".join(anc_names),
            "--pgen-prefixes", ",".join(pgen_prefixes),
            "--out-dirs", ",".join(out_dirs),
            "--corr-matrix", json.dumps(sim.get("corr_matrix", [])),
            "--heritability", str(sim.get("heritability", SIM_CFG.get("heritability", 0.4))),
            "--maf", str(sim.get("maf", SIM_CFG.get("maf", 0.05))),
            "--seed", str(sim.get("seed", SIM_CFG.get("seed", 42))),
            "--n_sims", str(sim.get("n_sims", SIM_CFG.get("n_sims", 10))),
            "--skip_thinning", str(sim.get("skip_thinning", SIM_CFG.get("skip_thinning", True))).lower(),
            "--thin_count_snps", str(sim.get("thin_count_snps", SIM_CFG.get("thin_count_snps", 1000000))),
            "--thin_count_inds", str(sim.get("thin_count_inds", SIM_CFG.get("thin_count_inds", 10000))),
        ], check=True)


rule computeSimGRM:
    container:
        "docker://gfanz/plink2:latest"
    conda:
        "../../envs/phenotypeSim.yml"
    threads: 4
    resources:
        nodes=1,
        mem_mb=16000,
        runtime=120,
    input:
        bed=OUT_DIR / "{ancestry}" / "simulations" / "{sim_name}" / "simulated.bed",
        bim=OUT_DIR / "{ancestry}" / "simulations" / "{sim_name}" / "simulated.bim",
        fam=OUT_DIR / "{ancestry}" / "simulations" / "{sim_name}" / "simulated.fam",
    output:
        grm_bin=OUT_DIR / "{ancestry}" / "simulations" / "{sim_name}" / "simulated.grm.bin",
        grm_id=OUT_DIR / "{ancestry}" / "simulations" / "{sim_name}" / "simulated.grm.id",
        grm_nbin=OUT_DIR / "{ancestry}" / "simulations" / "{sim_name}" / "simulated.grm.N.bin",
    params:
        prefix=lambda w, input: str(input.bed)[:-4],
    shell:
        """
        plink2 --bfile {params.prefix} --make-grm-bin --out {params.prefix}
        """


rule computeSimPCA:
    container:
        "docker://gfanz/plink2:latest"
    conda:
        "../../envs/phenotypeSim.yml"
    threads: 4
    resources:
        nodes=1,
        mem_mb=16000,
        runtime=120,
    input:
        bed=OUT_DIR / "{ancestry}" / "simulations" / "{sim_name}" / "simulated.bed",
        bim=OUT_DIR / "{ancestry}" / "simulations" / "{sim_name}" / "simulated.bim",
        fam=OUT_DIR / "{ancestry}" / "simulations" / "{sim_name}" / "simulated.fam",
    output:
        eigenvec=OUT_DIR / "{ancestry}" / "simulations" / "{sim_name}" / "simulated.eigenvec",
    params:
        prefix=lambda w, input: str(input.bed)[:-4],
    shell:
        """
        plink2 --bfile {params.prefix} --pca --out {params.prefix}
        """


rule extractSimPheno:
    container:
        "docker://gfanz/plink2:latest"
    conda:
        "../../envs/phenotypeSim.yml"
    threads: 1
    resources:
        nodes=1,
        mem_mb=4000,
        runtime=60,
    input:
        fam=OUT_DIR / "{ancestry}" / "simulations" / "{sim_name}" / "simulated.fam",
    output:
        pheno=OUT_DIR / "{ancestry}" / "simulations" / "{sim_name}" / "simulated_pheno1.pheno",
    shell:
        """
        awk 'BEGIN{{OFS=" "}}{{print $1, $2, $6}}' {input.fam} > {output.pheno}
        """