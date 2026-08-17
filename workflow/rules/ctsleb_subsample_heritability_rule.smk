SIM_CONFIG = config.get("phenotypeSimulation", {})
SNP_HERIT_CONFIG = config.get("snpHerit", {})

sim_ancestries = SIM_CONFIG.get("ancestries", ["AFR", "EUR"])
ANC1 = sim_ancestries[0]
ANC2 = sim_ancestries[1]
SIM_OUT_DIR = OUT_DIR / "simulations" / f"{ANC1}_{ANC2}"


def use_heritability_subsample():
    return bool(SNP_HERIT_CONFIG.get("subsample_before_heritability", False))


def heritability_sample_count():
    return int(SNP_HERIT_CONFIG.get("subsample_n", 20000))


def heritability_sample_seed():
    return int(SNP_HERIT_CONFIG.get("subsample_seed", 42))


def herit_prefix(anc):
    if use_heritability_subsample():
        return SIM_OUT_DIR / f"{anc}_simulation_herit"
    return SIM_OUT_DIR / f"{anc}_simulation"


rule subsampleSimForHeritability:
    log:
        OUT_DIR / "logs" / "subsampleSimForHeritability_{anc}.log",
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
    conda:
        "../../envs/ancNreport.yml"
    threads: 4
    resources:
        nodes=1,
        mem_mb=32000,
        runtime=240,
    input:
        bed=SIM_OUT_DIR / "{anc}_simulation.bed",
        bim=SIM_OUT_DIR / "{anc}_simulation.bim",
        fam=SIM_OUT_DIR / "{anc}_simulation.fam",
    output:
        bed=SIM_OUT_DIR / "{anc}_simulation_herit.bed",
        bim=SIM_OUT_DIR / "{anc}_simulation_herit.bim",
        fam=SIM_OUT_DIR / "{anc}_simulation_herit.fam",
    params:
        input_prefix=lambda wildcards, input: str(input.bed)[:-4],
        output_prefix=lambda wildcards, output: str(output.bed)[:-4],
        enabled=use_heritability_subsample(),
        n=heritability_sample_count(),
        seed=heritability_sample_seed(),
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname {output.bed})"

        if [[ "{params.enabled}" == "True" ]]; then
            echo "Subsampling {wildcards.anc} simulation to {params.n} individuals"
            plink2 --bfile {params.input_prefix} \
                --thin-indiv-count {params.n} \
                --seed {params.seed} \
                --make-bed \
                --out {params.output_prefix}
        else
            echo "subsample_before_heritability is false; copying full simulation for {wildcards.anc}"
            cp {input.bed} {output.bed}
            cp {input.bim} {output.bim}
            cp {input.fam} {output.fam}
        fi
        """


rule pruneSimVariants:
    log:
        OUT_DIR / "logs" / "pruneSimVariants_{anc}.log",
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
    conda:
        "../../envs/ancNreport.yml"
    threads: 8
    resources:
        nodes=1,
        mem_mb=32000,
        runtime=480,
    input:
        bed=rules.subsampleSimForHeritability.output.bed,
        bim=rules.subsampleSimForHeritability.output.bim,
        fam=rules.subsampleSimForHeritability.output.fam,
    output:
        bed=SIM_OUT_DIR / "{anc}_simulation_herit_pruned.bed",
        bim=SIM_OUT_DIR / "{anc}_simulation_herit_pruned.bim",
        fam=SIM_OUT_DIR / "{anc}_simulation_herit_pruned.fam",
    params:
        input_prefix=lambda wildcards, input: str(input.bed)[:-4],
        output_prefix=lambda wildcards, output: str(output.bed)[:-4],
        window_kb=int(SNP_HERIT_CONFIG.get("ld_window_kb", 500)),
        step=int(SNP_HERIT_CONFIG.get("ld_step", 10)),
        r2=float(SNP_HERIT_CONFIG.get("ld_r2", 0.1)),
    shell:
        """
        set -euo pipefail
        plink2 --bfile {params.input_prefix} \
            --indep-pairwise {params.window_kb} {params.step} {params.r2} \
            --out {params.output_prefix}

        plink2 --bfile {params.input_prefix} \
            --extract {params.output_prefix}.prune.in \
            --make-bed \
            --out {params.output_prefix}
        """


rule generateSimPCA:
    log:
        OUT_DIR / "logs" / "generateSimPCA_{anc}.log",
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
    threads: 8
    resources:
        nodes=1,
        mem_mb=64000,
        runtime=1440,
    input:
        bed=SIM_OUT_DIR / "{anc}_simulation_herit_pruned.bed",
        bim=SIM_OUT_DIR / "{anc}_simulation_herit_pruned.bim",
        fam=SIM_OUT_DIR / "{anc}_simulation_herit_pruned.fam",
    output:
        grm=SIM_OUT_DIR / "{anc}_simulation.grm.bin",
        grmid=SIM_OUT_DIR / "{anc}_simulation.grm.id",
        grmN=SIM_OUT_DIR / "{anc}_simulation.grm.N.bin",
        eigenvec=SIM_OUT_DIR / "{anc}_simulation.eigenvec",
    params:
        prefix=lambda wildcards, input: str(input.bed)[:-4],
        out_prefix=lambda wildcards: str(SIM_OUT_DIR / f"{wildcards.anc}_simulation"),
        npc=SNP_HERIT_CONFIG.get("npc", 10),
    shell:
        """
        set -euo pipefail
        plink2 --bfile {params.prefix} \
            --make-grm-bin \
            --pca approx {params.npc} \
            --out {params.out_prefix}
        """
