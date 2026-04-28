SIM_CONFIG = config.get("phenotypeSimulation", {})
sim_ancestries = SIM_CONFIG.get("ancestries", ["AFR", "EUR"])
ANC1 = sim_ancestries[0]
ANC2 = sim_ancestries[1]
SIM_OUT_DIR = OUT_DIR / "simulations" / f"{ANC1}_{ANC2}"


rule generateSimPCA:
    log:
        OUT_DIR / "logs" / "generateSimPCA_{anc}.log",
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/phenotypesim:latest"
    threads: 8
    resources:
        nodes=1,
        mem_mb=16000,
        runtime=60,
    input:
        bed=SIM_OUT_DIR / "{anc}_simulation.bed",
        bim=SIM_OUT_DIR / "{anc}_simulation.bim",
        fam=SIM_OUT_DIR / "{anc}_simulation.fam",
    output:
        eigenvec=SIM_OUT_DIR / "{anc}_simulation.eigenvec",
        grm_bin=SIM_OUT_DIR / "{anc}_simulation.grm.bin",
        grm_id=SIM_OUT_DIR / "{anc}_simulation.grm.id",
        grm_n=SIM_OUT_DIR / "{anc}_simulation.grm.N.bin",
    params:
        prefix=lambda wildcards, input: input.bed[:-4],
        npc=config.get("snpHerit", {}).get("npc", 10),
    shell:
        """
        plink2 --bfile {params.prefix} --make-grm-bin --pca approx {params.npc} --out {params.prefix}
        """


rule prepareSimPheno:
    container:
        "docker://alpine:latest"
    log:
        OUT_DIR / "logs" / "prepareSimPheno_{anc}.log",
    input:
        fam=SIM_OUT_DIR / "{anc}_simulation.fam",
    output:
        pheno=SIM_OUT_DIR / "{anc}_simulation_pheno1.pheno",
    shell:
        """
        awk '{{print $1, $2, $6}}' {input.fam} > {output.pheno}
        """


rule snpHeritSim:
    log:
        OUT_DIR / "logs" / "snpHeritSim_{anc}.log",
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/mash:latest"
    threads: 1
    resources:
        nodes=1,
        mem_mb=16000,
        runtime=480,
    input:
        grm=SIM_OUT_DIR / "{anc}_simulation.grm.bin",
        grmid=SIM_OUT_DIR / "{anc}_simulation.grm.id",
        grmN=SIM_OUT_DIR / "{anc}_simulation.grm.N.bin",
        eigen=SIM_OUT_DIR / "{anc}_simulation.eigenvec",
        pheno=SIM_OUT_DIR / "{anc}_simulation_pheno1.pheno",
    output:
        estimates=SIM_OUT_DIR / "{anc}_simulation_pheno1.estimates",
    params:
        method=config.get("snpHerit", {}).get("method", "AdjHE"),
        npc=config.get("snpHerit", {}).get("npc", 10),
        grm_prefix=lambda wildcards, input: input.grm[:-4],
        covar=config.get("snpHerit", {}).get("covar", ""),
    shell:
        """
        MASH --PC {input.eigen} \
          --prefix {params.grm_prefix} \
          --pheno {input.pheno} \
          --out {output.estimates} \
          --npc {params.npc} \
          --Method {params.method} \
          --covar {params.covar}
        """
