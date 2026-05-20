BUILD = config.get("ancestry", {}).get("build", "GRCh38")


if not INPUT_IS_PER_CHROMOSOME:
    rule crossmapStudyToB38:
        container:
            "oras://ghcr.io/coffm049/gdcgenomicsqc/crossmap:latest"
        conda:
            "../../envs/crossmap.yml"
        threads: 8
        resources:
            nodes=1,
            mem_mb=32000,
            runtime=120,
        input:
            pgen=OUT_DIR / "{subset}" / "f1.pgen",
            pvar=OUT_DIR / "{subset}" / "f1.pvar",
            psam=OUT_DIR / "{subset}" / "f1.psam",
            LDpvar=OUT_DIR / "{subset}" / "f1.ldpruned.pvar",
            chain=ancient(REF / "CrossMap" / "hg19ToHg38.over.chain.gz"),
        output:
            pgen=OUT_DIR / "{subset}" / "f1.b38.pgen",
            pvar=OUT_DIR / "{subset}" / "f1.b38.pvar",
            psam=OUT_DIR / "{subset}" / "f1.b38.psam",
            LDpgen=OUT_DIR / "{subset}" / "f1.b38.ldpruned.pgen",
            LDpvar=OUT_DIR / "{subset}" / "f1.b38.ldpruned.pvar",
            LDpsam=OUT_DIR / "{subset}" / "f1.b38.ldpruned.psam",
            tempDir=temp(directory(OUT_DIR / "{subset}" / "intermediates" / "crossmap")),
        params:
            input_prefix=OUT_DIR / "{subset}" / "f1",
            output_prefix=OUT_DIR / "{subset}" / "f1.b38",
            chain=ancient(REF / "CrossMap" / "hg19ToHg38.over.chain.gz"),
        run:
            if BUILD == "GRCh38":
                shell("""
                    mkdir -p {output.tempDir}
                    plink2 --pfile {params.input_prefix} \
                           --make-pgen \
                           --out {params.output_prefix}
                """)
            else:
                shell("""
                    mkdir -p {output.tempDir}

                    plink2 --pfile {params.input_prefix} --make-bed --out {output.tempDir}/study

                    awk '{{print $1, $4-1, $4, $2}}' {output.tempDir}/study.bim > {output.tempDir}/study_pos.bed

                    CrossMap.py bed {params.chain} {output.tempDir}/study_pos.bed {output.tempDir}/study_hg38

                    awk '{{print $4}}' {output.tempDir}/study_hg38.bed > {output.tempDir}/lifted_snps.txt
                    awk '{{print $4, $3}}' {output.tempDir}/study_hg38.bed > {output.tempDir}/new_pos.txt
                    awk '{{print $4, $1}}' {output.tempDir}/study_hg38.bed > {output.tempDir}/new_chr.txt

                    plink2 --bfile {output.tempDir}/study --extract {output.tempDir}/lifted_snps.txt --make-bed --out {output.tempDir}/step1
                    plink2 --bfile {output.tempDir}/step1 --update-map {output.tempDir}/new_pos.txt --make-bed --out {output.tempDir}/step2
                    plink2 --bfile {output.tempDir}/step2 --update-chr {output.tempDir}/new_chr.txt --make-pgen --out {params.output_prefix}
                """)
            shell("""
                awk 'NR>1 {{print $2}}' {params.input_prefix}.ldpruned.pvar > {output.tempDir}/ldpruned_vars.txt
                plink2 --pfile {params.output_prefix} --extract {output.tempDir}/ldpruned_vars.txt --make-pgen --out {params.output_prefix}.ldpruned
            """)
