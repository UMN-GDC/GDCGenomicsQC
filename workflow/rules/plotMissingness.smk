rule plotMissingness:
    log:
        OUT_DIR / "logs" / "plotMissingness_{subset}.log",
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
    conda:
        "../../envs/ancNreport.yml"
    threads: 1
    resources:
        nodes=1,
        mem_mb=8000,
        runtime=30,
    output:
        smissIMG=report(
            OUT_DIR / "{subset}" / "figures" / "smiss.svg",
            caption="../../report/smiss.rst",
            category="Quality Control",
        ),
        vmissIMG=report(
            OUT_DIR / "{subset}" / "figures" / "vmiss.svg",
            caption="../../report/vmiss.rst",
            category="Quality Control",
        ),
    input:
        smiss=expand(OUT_DIR / "{{subset}}" / "initial_{CHR}.smiss", CHR=CHROMOSOMES),
        vmiss=expand(OUT_DIR / "{{subset}}" / "initial_{CHR}.vmiss", CHR=CHROMOSOMES),
    params:
        smiss_files=lambda wildcards, input: " ".join(input.smiss),
        vmiss_files=lambda wildcards, input: " ".join(input.vmiss),
    shell:
        """
    Rscript scripts/plotMissingness.R "{params.smiss_files}" "{params.vmiss_files}" {output.smissIMG} {output.vmissIMG}
    """
