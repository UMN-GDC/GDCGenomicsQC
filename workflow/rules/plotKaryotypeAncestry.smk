def get_local_ancestry_samples():
    sis_file = OUT_DIR / "02-localAncestry" / "chr20.lai.sis.tsv"
    if not sis_file.exists():
        return []
    return pd.read_csv(sis_file, header=None)[0].tolist()


SAMPLES = get_local_ancestry_samples()


checkpoint generateKaryotypeAncestryPlots:
    log:
        OUT_DIR / "logs" / "generateKaryotypeAncestryPlots.log",
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
    conda:
        "../../envs/karyoploteR.yml"
    resources:
        nodes=1,
        mem_mb=16000,
        runtime=60,
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
    params:
        msp_dir=OUT_DIR / "02-localAncestry",
        figures_dir=lambda wildcards: OUT_DIR
        / "02-localAncestry"
        / config.get("localAncestry", {}).get("figures", "figures"),
        chromosomes="1-22",
    shell:
        """
        mkdir -p {params.figures_dir}
        for SAMPLE in {" ".join(SAMPLES)}; do
            Rscript ../../scripts/plotKaryotypeAncestry.R \
                --msp-dir {params.msp_dir} \
                --sample $SAMPLE \
                --chromosomes {params.chromosomes} \
                --output {params.figures_dir}/$SAMPLE'_karyotype.pdf'
        done
        """
