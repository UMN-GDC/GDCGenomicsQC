SIM_CONFIG = config.get("phenotypeSimulation", {})

sim_ancestries = SIM_CONFIG.get("ancestries", ["AFR", "EUR"])
if len(sim_ancestries) != 2:
    raise ValueError("phenotypeSimulation.ancestries must be exactly 2 values")

ANC1 = sim_ancestries[0]
ANC2 = sim_ancestries[1]

SIM_OUT_DIR = OUT_DIR / "simulations" / f"{ANC1}_{ANC2}"


rule simulatePhenotype:
    log:
        OUT_DIR / "logs" / "simulatePhenotype.log",
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/phenotypesim:latest"
    conda:
        "../../envs/phenotypeSim.yml"
    threads: 8
    resources:
        nodes=1,
        mem_mb=32000,
        runtime=240,
    input:
        anc1_pgen=OUT_DIR / ANC1 / "initialFilter.pgen",
        anc1_pvar=OUT_DIR / ANC1 / "initialFilter.pvar",
        anc1_psam=OUT_DIR / ANC1 / "initialFilter.psam",
        anc2_pgen=OUT_DIR / ANC2 / "initialFilter.pgen",
        anc2_pvar=OUT_DIR / ANC2 / "initialFilter.pvar",
        anc2_psam=OUT_DIR / ANC2 / "initialFilter.psam",
    output:
        sim_dir=directory(SIM_OUT_DIR),
        anc1_fam=SIM_OUT_DIR / f"{ANC1}_simulation.fam",
        anc2_fam=SIM_OUT_DIR / f"{ANC2}_simulation.fam",
        anc1_bed=SIM_OUT_DIR / f"{ANC1}_simulation.bed",
        anc1_bim=SIM_OUT_DIR / f"{ANC1}_simulation.bim",
        anc2_bed=SIM_OUT_DIR / f"{ANC2}_simulation.bed",
        anc2_bim=SIM_OUT_DIR / f"{ANC2}_simulation.bim",
    params:
        n_sims=SIM_CONFIG.get("n_sims", 10),
        heritability=SIM_CONFIG.get("heritability", 0.4),
        rho=SIM_CONFIG.get("rho", 0.8),
        maf=SIM_CONFIG.get("maf", 0.05),
        seed=SIM_CONFIG.get("seed", 42),
        skip_thinning=SIM_CONFIG.get("skip_thinning", True),
        thin_count_snps=SIM_CONFIG.get("thin_count_snps", 1000000),
        thin_count_inds=SIM_CONFIG.get("thin_count_inds", 10000),
        anc1=ANC1,
        anc2=ANC2,
        script_dir=Path(workflow.basedir) / "scripts",
    shell:
        """
        set -euo pipefail

        anc1_pfile="{input.anc1_pgen}"
        anc1_pfile="${{anc1_pfile%.pgen}}"
        anc2_pfile="{input.anc2_pgen}"
        anc2_pfile="${{anc2_pfile%.pgen}}"

        Rscript {params.script_dir}/runPhenotypeSimulation.R \
            --ancestry1 "$anc1_pfile" \
            --ancestry2 "$anc2_pfile" \
            --out_dir {output.sim_dir} \
            --anc1_name {params.anc1} \
            --anc2_name {params.anc2} \
            --n_sims {params.n_sims} \
            --seed {params.seed} \
            --heritability {params.heritability} \
            --rho {params.rho} \
            --maf {params.maf} \
            --skip_thinning {params.skip_thinning} \
            --thin_count_snps {params.thin_count_snps} \
            --thin_count_inds {params.thin_count_inds}
        """
