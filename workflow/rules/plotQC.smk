if INPUT_IS_PER_CHROMOSOME:
    rule plotHWE:
        log:
            OUT_DIR / "logs" / "plotHWE_{subset}_{CHR}.log",
        container:
            "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
        conda:
            "../../envs/ancNreport.yml"
        envmodules: *[m for m in (config.get("plink_module"), config.get("R_module")) if m]
        threads: 4
        resources:
            nodes=1,
            mem_mb=16000,
            runtime=30,
        input:
            hardy=OUT_DIR / "{subset}" / "f1.f2_{CHR}.hardy",
        output:
            plot=report(
                OUT_DIR / "{subset}" / "figures" / "hwe_histogram_{CHR}.svg",
                caption="Histogram of -log10 Hardy-Weinberg equilibrium p-values",
                category="Quality Control",
            ),
        params:
            scripts_dir=SCRIPTS_DIR,
        shell:
            """
            mkdir -p "$(dirname {output.plot})"
            if [ -f {input.hardy} ]; then
                Rscript {params.scripts_dir}/plotHWE.R {input.hardy} {output.plot}
            else
                echo "Warning: {input.hardy} not found, skipping HWE plot" >> {log}
            fi
            """

    rule plotHeterozygosity:
        log:
            OUT_DIR / "logs" / "plotHeterozygosity_{subset}_{CHR}.log",
        container:
            "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
        conda:
            "../../envs/ancNreport.yml"
        envmodules: *[m for m in (config.get("plink_module"), config.get("R_module")) if m]
        threads: 4
        resources:
            nodes=1,
            mem_mb=16000,
            runtime=30,
        input:
            het=OUT_DIR / "{subset}" / "f1.f2_{CHR}.het",
        output:
            plot=report(
                OUT_DIR / "{subset}" / "figures" / "het_distribution_{CHR}.svg",
                caption="Distribution of inbreeding coefficients (F) after standard QC",
                category="Quality Control",
            ),
        params:
            scripts_dir=SCRIPTS_DIR,
        shell:
            """
            mkdir -p "$(dirname {output.plot})"
            if [ -f {input.het} ]; then
                Rscript {params.scripts_dir}/plotHeterozygosity.R {input.het} {output.plot}
            else
                echo "Warning: {input.het} not found, skipping heterozygosity plot" >> {log}
            fi
            """

else:
    rule plotHWE:
        log:
            OUT_DIR / "logs" / "plotHWE_{subset}.log",
        container:
            "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
        conda:
            "../../envs/ancNreport.yml"
        envmodules: *[m for m in (config.get("plink_module"), config.get("R_module")) if m]
        threads: 4
        resources:
            nodes=1,
            mem_mb=16000,
            runtime=30,
        input:
            hardy=OUT_DIR / "{subset}" / "f1.b38.f2.hardy",
        output:
            plot=report(
                OUT_DIR / "{subset}" / "figures" / "hwe_histogram.svg",
                caption="Histogram of -log10 Hardy-Weinberg equilibrium p-values",
                category="Quality Control",
            ),
        params:
            scripts_dir=SCRIPTS_DIR,
        shell:
            """
            mkdir -p "$(dirname {output.plot})"
            if [ -f {input.hardy} ]; then
                Rscript {params.scripts_dir}/plotHWE.R {input.hardy} {output.plot}
            else
                echo "Warning: {input.hardy} not found, skipping HWE plot" >> {log}
            fi
            """

    rule plotHeterozygosity:
        log:
            OUT_DIR / "logs" / "plotHeterozygosity_{subset}.log",
        container:
            "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
        conda:
            "../../envs/ancNreport.yml"
        envmodules: *[m for m in (config.get("plink_module"), config.get("R_module")) if m]
        threads: 4
        resources:
            nodes=1,
            mem_mb=16000,
            runtime=30,
        input:
            het=OUT_DIR / "{subset}" / "f1.b38.f2.het",
        output:
            plot=report(
                OUT_DIR / "{subset}" / "figures" / "het_distribution.svg",
                caption="Distribution of inbreeding coefficients (F) after standard QC",
                category="Quality Control",
            ),
        params:
            scripts_dir=SCRIPTS_DIR,
        shell:
            """
            mkdir -p "$(dirname {output.plot})"
            if [ -f {input.het} ]; then
                Rscript {params.scripts_dir}/plotHeterozygosity.R {input.het} {output.plot}
            else
                echo "Warning: {input.het} not found, skipping heterozygosity plot" >> {log}
            fi
            """

    rule plotRelatedness:
        log:
            OUT_DIR / "logs" / "plotRelatedness_{subset}.log",
        container:
            "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
        conda:
            "../../envs/ancNreport.yml"
        envmodules: *[m for m in (config.get("plink_module"), config.get("R_module")) if m]
        threads: 8
        resources:
            nodes=1,
            mem_mb=32000,
            runtime=60,
        input:
            king=OUT_DIR / "{subset}" / "f1.b38.ldpruned.unrelated_grm.king",
        output:
            plot=report(
                OUT_DIR / "{subset}" / "figures" / "relatedness_histogram.svg",
                caption="Histogram of pairwise KING kinship coefficients (unrelated set)",
                category="Quality Control",
            ),
        params:
            scripts_dir=SCRIPTS_DIR,
        shell:
            """
            mkdir -p "$(dirname {output.plot})"
            if [ -f {input.king} ]; then
                Rscript {params.scripts_dir}/plotRelatedness.R {input.king} {output.plot}
            else
                echo "Warning: {input.king} not found, skipping relatedness plot" >> {log}
            fi
            """
