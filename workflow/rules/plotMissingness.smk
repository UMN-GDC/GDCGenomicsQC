if INPUT_IS_PER_CHROMOSOME:
    rule plotSampleVariantMissingness:
        log:
            OUT_DIR / "logs" / "plotSampleVariantMissingness_{subset}_{CHR}.log",
        container:
            "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
        conda:
            "../../envs/ancNreport.yml"
        envmodules: *([config.get("R_module")] if config.get("R_module") else [])
        threads: 1
        resources:
            nodes=1,
            mem_mb=8000,
            runtime=30,
        output:
            smissIMG=report(
                OUT_DIR / "{subset}" / "figures" / "smiss_{CHR}.svg",
                caption="../../report/smiss.rst",
                category="Quality Control",
            ),
            vmissIMG=report(
                OUT_DIR / "{subset}" / "figures" / "vmiss_{CHR}.svg",
                caption="../../report/vmiss.rst",
                category="Quality Control",
            ),
        input:
            smiss=OUT_DIR / "{subset}" / "initial_{CHR}.smiss",
            vmiss=OUT_DIR / "{subset}" / "initial_{CHR}.vmiss",
        params:
            scripts_dir=SCRIPTS_DIR,
        shell:
            """
        Rscript {params.scripts_dir}/plotMissingness.R "{input.smiss}" "{input.vmiss}" {output.smissIMG} {output.vmissIMG}
        """

else:
    rule plotSampleVariantMissingness:
        log:
            OUT_DIR / "logs" / "plotSampleVariantMissingness_{subset}.log",
        container:
            "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
        conda:
            "../../envs/ancNreport.yml"
        envmodules: *([config.get("R_module")] if config.get("R_module") else [])
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
            smiss=OUT_DIR / "{subset}" / "initial.smiss",
            vmiss=OUT_DIR / "{subset}" / "initial.vmiss",
        params:
            scripts_dir=SCRIPTS_DIR,
        shell:
            """
        Rscript {params.scripts_dir}/plotMissingness.R "{input.smiss}" "{input.vmiss}" {output.smissIMG} {output.vmissIMG}
        """
