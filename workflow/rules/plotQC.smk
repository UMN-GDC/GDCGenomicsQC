if not INPUT_IS_PER_CHROMOSOME:

    rule plotHWE:
        log:
            OUT_DIR / "logs" / "plotHWE_{subset}.log",
        container:
            "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
        conda:
            "../../envs/ancNreport.yml"
        envmodules:
            mod("plink2", "R")
        threads: 4
        resources:
            nodes=1,
            mem_mb=16000,
            runtime=30,
        input:
            pgen=OUT_DIR / "{subset}" / "f1.b38.f2.pgen",
            pvar=OUT_DIR / "{subset}" / "f1.b38.f2.pvar",
            psam=OUT_DIR / "{subset}" / "f1.b38.f2.psam",
        output:
            plot=report(
                OUT_DIR / "{subset}" / "figures" / "hwe_histogram.svg",
                caption="Histogram of -log10 Hardy-Weinberg equilibrium p-values",
                category="Quality Control",
            ),
        params:
            prefix=lambda wildcards, input: str(input.pgen)[:-5],
            scripts_dir=SCRIPTS_DIR,
        shell:
            """
            mkdir -p "$(dirname {output.plot})"
            plink2 --pfile {params.prefix} --hardy --out {params.prefix}_hwe --threads {threads}
            Rscript {params.scripts_dir}/plotHWE.R {params.prefix}_hwe.hwe {output.plot}
            """


    rule plotHeterozygosity:
        log:
            OUT_DIR / "logs" / "plotHeterozygosity_{subset}.log",
        container:
            "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
        conda:
            "../../envs/ancNreport.yml"
        envmodules:
            mod("plink2", "R")
        threads: 4
        resources:
            nodes=1,
            mem_mb=16000,
            runtime=30,
        input:
            pgen=OUT_DIR / "{subset}" / "f1.b38.f2.pgen",
            pvar=OUT_DIR / "{subset}" / "f1.b38.f2.pvar",
            psam=OUT_DIR / "{subset}" / "f1.b38.f2.psam",
        output:
            plot=report(
                OUT_DIR / "{subset}" / "figures" / "het_distribution.svg",
                caption="Distribution of inbreeding coefficients (F) after standard QC",
                category="Quality Control",
            ),
        params:
            prefix=lambda wildcards, input: str(input.pgen)[:-5],
            scripts_dir=SCRIPTS_DIR,
        shell:
            """
            mkdir -p "$(dirname {output.plot})"
            mkdir -p {params.prefix}_het_tmp
            plink2 --pfile {params.prefix} \
                --indep-pairwise 50 5 0.2 \
                --out {params.prefix}_het_tmp/indep \
                --threads {threads}
            plink2 --pfile {params.prefix} \
                --extract {params.prefix}_het_tmp/indep.prune.in \
                --het \
                --out {params.prefix}_het \
                --threads {threads}
            Rscript {params.scripts_dir}/plotHeterozygosity.R {params.prefix}_het.het {output.plot}
            """


    rule plotRelatedness:
        log:
            OUT_DIR / "logs" / "plotRelatedness_{subset}.log",
        container:
            "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
        conda:
            "../../envs/ancNreport.yml"
        envmodules:
            mod("plink2", "R")
        threads: 8
        resources:
            nodes=1,
            mem_mb=32000,
            runtime=60,
        input:
            pgen=OUT_DIR / "{subset}" / "f1.b38.ldpruned.unrelated.pgen",
            pvar=OUT_DIR / "{subset}" / "f1.b38.ldpruned.unrelated.pvar",
            psam=OUT_DIR / "{subset}" / "f1.b38.ldpruned.unrelated.psam",
        output:
            plot=report(
                OUT_DIR / "{subset}" / "figures" / "relatedness_histogram.svg",
                caption="Histogram of pairwise KING kinship coefficients (unrelated set)",
                category="Quality Control",
            ),
        params:
            prefix=lambda wildcards, input: str(input.pgen)[:-5],
            scripts_dir=SCRIPTS_DIR,
        shell:
            """
            mkdir -p "$(dirname {output.plot})"
            plink2 --pfile {params.prefix} \
                --make-king \
                --out {params.prefix}_kin \
                --threads {threads}
            Rscript {params.scripts_dir}/plotRelatedness.R {params.prefix}_kin.kin0 {output.plot}
            """
