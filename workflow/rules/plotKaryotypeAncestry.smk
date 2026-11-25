rule generateKaryotypeAncestryPlots:
    input:
        expand(OUT_DIR / "02-localAncestry" / "chr{CHR}.lai.msp.tsv", CHR=CHROMOSOMES),
    output:
        directory(
            OUT_DIR
            / "02-localAncestry"
            / config.get("localAncestry", {}).get("figures", "figures")
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
        scripts_dir=SCRIPTS_DIR,
    shell:
        mkdir -p {output}
        # Discover samples from existing .lai.sis.tsv files
        SAMPLES=""
        for sis_file in {params.msp_dir}/chr*.lai.sis.tsv; do
            [ -f "$sis_file" ] || continue
            SAMPLE=$(basename "$sis_file" .lai.sis.tsv)
            SAMPLES="$SAMPLES $SAMPLE"
        done
        for SAMPLE in $SAMPLES; do
            [ -z "$SAMPLE" ] && continue
            Rscript {params.scripts_dir}/plotKaryotypeAncestry.R \
                --msp-dir {params.msp_dir} \
                --sample $SAMPLE \
                --chromosomes {params.chromosomes} \
                --output {output}/$SAMPLE'_karyotype.pdf'
        done