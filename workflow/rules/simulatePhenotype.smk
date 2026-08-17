from pathlib import Path

SIM_CONFIG = config.get("phenotypeSimulation", {})

sim_ancestries = SIM_CONFIG.get("ancestries", ["AFR", "EUR"])
if len(sim_ancestries) != 2:
    raise ValueError("phenotypeSimulation.ancestries must be exactly 2 values")

ANC1 = sim_ancestries[0]
ANC2 = sim_ancestries[1]
SIM_OUT_DIR = OUT_DIR / "simulations" / f"{ANC1}_{ANC2}"
SIM_INPUT_PREFIXES = SIM_CONFIG.get("input_prefixes", {})


def get_sim_input_prefix(anc):
    return str(SIM_INPUT_PREFIXES.get(anc, OUT_DIR / anc / "initialFilter"))


def get_sim_input_kind(anc):
    prefix = get_sim_input_prefix(anc)

    pfile_exts = ("pgen", "pvar", "psam")
    if all(Path(f"{prefix}.{ext}").exists() for ext in pfile_exts):
        return "pfile"

    bfile_exts = ("bed", "bim", "fam")
    if all(Path(f"{prefix}.{ext}").exists() for ext in bfile_exts):
        return "bfile"

    if anc in SIM_INPUT_PREFIXES:
        raise FileNotFoundError(
            f"No supported PLINK input set found for ancestry {anc!r} at prefix {prefix}. "
            "Expected either .pgen/.pvar/.psam or .bed/.bim/.fam."
        )

    return "pfile"


def get_sim_input_file(anc, ext):
    prefix = get_sim_input_prefix(anc)
    kind = get_sim_input_kind(anc)

    if kind == "pfile" and ext in {"pgen", "pvar", "psam"}:
        return f"{prefix}.{ext}"
    if kind == "bfile" and ext in {"bed", "bim", "fam"}:
        return f"{prefix}.{ext}"

    return []


rule simulateBivariatePhenotypes:
    log:
        OUT_DIR / "logs" / "simulateBivariatePhenotypes.log",
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
        anc1_pgen=lambda wildcards: get_sim_input_file(ANC1, "pgen"),
        anc1_pvar=lambda wildcards: get_sim_input_file(ANC1, "pvar"),
        anc1_psam=lambda wildcards: get_sim_input_file(ANC1, "psam"),
        anc1_bed=lambda wildcards: get_sim_input_file(ANC1, "bed"),
        anc1_bim=lambda wildcards: get_sim_input_file(ANC1, "bim"),
        anc1_fam=lambda wildcards: get_sim_input_file(ANC1, "fam"),
        anc2_pgen=lambda wildcards: get_sim_input_file(ANC2, "pgen"),
        anc2_pvar=lambda wildcards: get_sim_input_file(ANC2, "pvar"),
        anc2_psam=lambda wildcards: get_sim_input_file(ANC2, "psam"),
        anc2_bed=lambda wildcards: get_sim_input_file(ANC2, "bed"),
        anc2_bim=lambda wildcards: get_sim_input_file(ANC2, "bim"),
        anc2_fam=lambda wildcards: get_sim_input_file(ANC2, "fam"),
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
        anc1_prefix=get_sim_input_prefix(ANC1),
        anc2_prefix=get_sim_input_prefix(ANC2),
        anc1_kind=get_sim_input_kind(ANC1),
        anc2_kind=get_sim_input_kind(ANC2),
        script_dir=Path(workflow.basedir) / "scripts",
    shell:
        """
        set -euo pipefail

        mkdir -p {output.sim_dir}

        anc1_pfile="{params.anc1_prefix}"
        anc2_pfile="{params.anc2_prefix}"

        if [[ "{params.anc1_kind}" == "bfile" ]]; then
            anc1_pfile="{output.sim_dir}/{params.anc1}_input"
            plink2 --bfile "/scratch.global/saonli/GDCGenomicsQC/CTSLEB/AFR/CTSLEB_AFR" \
  	--set-all-var-ids 'chr@:#:$r:$a' \
  	--rm-dup force-first \
  	--make-pgen \
  	--out "$anc1_pfile"
        fi

        if [[ "{params.anc2_kind}" == "bfile" ]]; then
            anc2_pfile="{output.sim_dir}/{params.anc2}_input"
         plink2 --bfile "/scratch.global/saonli/GDCGenomicsQC/CTSLEB/EUR/CTSLEB_EUR" \
  	--set-all-var-ids 'chr@:#:$r:$a' \
  	--rm-dup force-first \
  	--make-pgen \
  	--out "$anc2_pfile"
        fi

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
