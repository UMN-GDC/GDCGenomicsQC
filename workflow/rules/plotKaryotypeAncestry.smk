def get_local_ancestry_samples():
    sis_file = OUT_DIR / "02-localAncestry" / "chr20.lai.sis.tsv"
    if not sis_file.exists():
        return []
    return pd.read_csv(sis_file, header=None)[0].tolist()


SAMPLES = get_local_ancestry_samples()
SAMPLES_STR = " ".join(SAMPLES)


checkpoint generateKaryotypeAncestryPlots:
    input:
        expand(OUT_DIR / "02-localAncestry" / "chr{CHR}.lai.msp.tsv", CHR=CHROMOSOMES),
    output:
        expand(
            OUT_DIR
            / "02-localAncestry"
            / config.get("localAncestry", {}).get("figures", "figures")
            / "{sample}_karyotype.pdf",
            sample=SAMPLES,
        ),
    log:
        OUT_DIR / "logs" / "generateKaryotypeAncestryPlots.log",
    conda:
        "../../envs/karyoploteR.yml"
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
    envmodules:
        *([config.get("R_module")] if config.get("R_module") else []),
    resources:
        nodes=1,
        mem_mb=16000,
        runtime=60,
    params:
        msp_dir=OUT_DIR / "02-localAncestry",
        figures_dir=OUT_DIR
        / "02-localAncestry"
        / config.get("localAncestry", {}).get("figures", "figures"),
        chromosomes="1-22",
        samples_str=SAMPLES_STR,
        scripts_dir=SCRIPTS_DIR,
    shell:
        """
        mkdir -p {params.figures_dir}
        for SAMPLE in {params.samples_str}; do
            Rscript {params.scripts_dir}/plotKaryotypeAncestry.R \
                --msp-dir {params.msp_dir} \
                --sample $SAMPLE \
                --chromosomes {params.chromosomes} \
                --output {params.figures_dir}/$SAMPLE'_karyotype.pdf'
        done
        """
