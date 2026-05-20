BUILD = config.get("ancestry", {}).get("build", "GRCh38")


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
        pgen=OUT_DIR / "full" / "f1.pgen",
        pvar=OUT_DIR / "full" / "f1.pvar",
        psam=OUT_DIR / "full" / "f1.psam",
        chain=ancient(REF / "CrossMap" / "hg19ToHg38.over.chain.gz"),
    output:
        pgen=OUT_DIR / "full" / "f1.b38.pgen",
        pvar=OUT_DIR / "full" / "f1.b38.pvar",
        psam=OUT_DIR / "full" / "f1.b38.psam",
        tempDir=temp(directory(OUT_DIR / "full" / "intermediates" / "crossmap")),
    params:
        input_prefix=OUT_DIR / "full" / "f1",
        output_prefix=OUT_DIR / "full" / "f1.b38",
        chain=ancient(REF / "CrossMap" / "hg19ToHg38.over.chain.gz"),
    run:
        if BUILD == "GRCh38":
            shell("""
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
