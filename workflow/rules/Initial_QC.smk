if INPUT_IS_PER_CHROMOSOME:
    rule concatLdPruned:
        log:
            OUT_DIR / "logs" / "concatLdPruned_{subset}.log",
        container:
            "docker://gfanz/plink2:latest"
        conda:
            "../../envs/ancNreport.yml"
        envmodules: *([config.get("plink_module")] if config.get("plink_module") else [])
        threads: 4
        resources:
            nodes=1,
            mem_mb=16000,
            runtime=60,
        output:
            LDpgen=OUT_DIR / "{subset}" / "f1.ldpruned.pgen",
            LDpvar=OUT_DIR / "{subset}" / "f1.ldpruned.pvar",
            LDpsam=OUT_DIR / "{subset}" / "f1.ldpruned.psam",
            tempDir=temp(
                directory(OUT_DIR / "{subset}" / "intermediates" / "ldpruned_concat")
            ),
        input:
            LDpgen=expand(
                OUT_DIR / "{{subset}}" / "f1.ldpruned_{CHR}.pgen", CHR=CHROMOSOMES
            ),
            LDpvar=expand(
                OUT_DIR / "{{subset}}" / "f1.ldpruned_{CHR}.pvar", CHR=CHROMOSOMES
            ),
            LDpsam=expand(
                OUT_DIR / "{{subset}}" / "f1.ldpruned_{CHR}.psam", CHR=CHROMOSOMES
            ),
        params:
            output_prefix=lambda wildcards, output: output.LDpgen[:-5],
        shell:
            """
            mkdir -p {output.tempDir}
            > {output.tempDir}/mergelist.txt
            for f in {input.LDpgen}; do
                echo "${{f%.pgen}}" >> {output.tempDir}/mergelist.txt
            done
            plink2 --pmerge-list {output.tempDir}/mergelist.txt \
                   --make-pgen \
                   --threads {threads} \
                   --out {params.output_prefix}
            """
