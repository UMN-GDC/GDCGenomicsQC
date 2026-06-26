if INPUT_IS_PER_CHROMOSOME:
    rule crossmapFullToB38:
        log:
            OUT_DIR / "logs" / "crossmapFullToB38_{subset}_{CHR}.log",
        container:
            "oras://ghcr.io/coffm049/gdcgenomicsqc/crossmap:latest"
        conda:
            "../../envs/crossmap.yml"
        envmodules: *([config.get("plink_module")] if config.get("plink_module") else [])
        threads: 2
        resources:
            nodes=1,
            mem_mb=8000,
            runtime=30,
        input:
            pgen=OUT_DIR / "{subset}" / "f1_{CHR}.pgen",
            pvar=OUT_DIR / "{subset}" / "f1_{CHR}.pvar",
            psam=OUT_DIR / "{subset}" / "f1_{CHR}.psam",
            chain=ancient(REF / "CrossMap" / "hg19ToHg38.over.chain.gz"),
        output:
            pgen=OUT_DIR / "{subset}" / "f1.b38_{CHR}.pgen",
            pvar=OUT_DIR / "{subset}" / "f1.b38_{CHR}.pvar",
            psam=OUT_DIR / "{subset}" / "f1.b38_{CHR}.psam",
            tempDir=temp(directory(OUT_DIR / "{subset}" / "intermediates" / "crossmap_{CHR}")),
        params:
            input_prefix=lambda wildcards: str(OUT_DIR / wildcards.subset / f"f1_{wildcards.CHR}"),
            output_prefix=lambda wildcards: str(OUT_DIR / wildcards.subset / f"f1.b38_{wildcards.CHR}"),
        run:
            _build = config.get("build", "GRCh38")
            if _build == "GRCh38":
                shell("""
                    mkdir -p {output.tempDir}
                    plink2 --pfile {params.input_prefix} \
                           --set-all-var-ids 'chr@:#:$r:$a' \
                           --make-pgen \
                           --out {params.output_prefix}
                """)
            else:
                shell("""
                    mkdir -p {output.tempDir}

                    awk '$1 !~ /^#/ {{print "chr"$1, $2-1, $2, $3}}' {params.input_prefix}.pvar > {output.tempDir}/study_pos.bed
                    N_INPUT=$(wc -l < {output.tempDir}/study_pos.bed)
                    CrossMap bed {params.chain} {output.tempDir}/study_pos.bed {output.tempDir}/study_hg38.bed
                    N_LIFTED=$(awk 'END{{print NR}}' {output.tempDir}/study_hg38.bed 2>/dev/null || echo 0)
                    N_DROPPED=$((N_INPUT - N_LIFTED))
                    echo "[CrossMap per-chr] Lifted $N_LIFTED / $N_INPUT variants ($N_DROPPED dropped)" >> {log} 2>&1
                    awk '{{print $4}}' {output.tempDir}/study_hg38.bed > {output.tempDir}/lifted_snps.txt
                    awk '{{print $4, $3}}' {output.tempDir}/study_hg38.bed > {output.tempDir}/new_pos.txt
                    awk '{{gsub(/^chr/,"",$1); print $4, $1}}' {output.tempDir}/study_hg38.bed > {output.tempDir}/new_chr.txt

                    plink2 --pfile {params.input_prefix} --extract {output.tempDir}/lifted_snps.txt --make-pgen --out {output.tempDir}/step1
                    plink2 --pfile {output.tempDir}/step1 --update-map {output.tempDir}/new_pos.txt --make-pgen --out {output.tempDir}/step2
                    plink2 --pfile {output.tempDir}/step2 --update-chr {output.tempDir}/new_chr.txt --sort-vars --set-all-var-ids 'chr@:#:$r:$a' --make-pgen --out {params.output_prefix}
                """)

else:
    rule crossmapFullToB38:
        log:
            OUT_DIR / "logs" / "crossmapFullToB38_{subset}.log",
        container:
            "oras://ghcr.io/coffm049/gdcgenomicsqc/crossmap:latest"
        conda:
            "../../envs/crossmap.yml"
        envmodules: *[m for m in (config.get("plink_module"), config.get("crossmap_module")) if m]
        threads: 8
        resources:
            nodes=1,
            mem_mb=32000,
            runtime=120,
        input:
            pgen=OUT_DIR / "{subset}" / "f1.pgen",
            pvar=OUT_DIR / "{subset}" / "f1.pvar",
            psam=OUT_DIR / "{subset}" / "f1.psam",
            chain=ancient(REF / "CrossMap" / "hg19ToHg38.over.chain.gz"),
        output:
            pgen=OUT_DIR / "{subset}" / "f1.b38.pgen",
            pvar=OUT_DIR / "{subset}" / "f1.b38.pvar",
            psam=OUT_DIR / "{subset}" / "f1.b38.psam",
            tempDir=temp(directory(OUT_DIR / "{subset}" / "intermediates" / "crossmap")),
        params:
            input_prefix=lambda wildcards: str(OUT_DIR / wildcards.subset / "f1"),
            output_prefix=lambda wildcards: str(OUT_DIR / wildcards.subset / "f1.b38"),
            chain=ancient(REF / "CrossMap" / "hg19ToHg38.over.chain.gz"),
        run:
            _build = config.get("build", "GRCh38")
            if _build == "GRCh38":
                shell("""
                    mkdir -p {output.tempDir}
                    plink2 --pfile {params.input_prefix} \
                           --set-all-var-ids 'chr@:#:$r:$a' \
                           --make-pgen \
                           --out {params.output_prefix}
                """)
            else:
                shell("""
                    mkdir -p {output.tempDir}

                    awk '$1 !~ /^#/ {{print "chr"$1, $2-1, $2, $3}}' {params.input_prefix}.pvar > {output.tempDir}/study_pos.bed
                    N_INPUT=$(wc -l < {output.tempDir}/study_pos.bed)
                    CrossMap bed {params.chain} {output.tempDir}/study_pos.bed {output.tempDir}/study_hg38.bed
                    N_LIFTED=$(awk 'END{{print NR}}' {output.tempDir}/study_hg38.bed 2>/dev/null || echo 0)
                    N_DROPPED=$((N_INPUT - N_LIFTED))
                    echo "[CrossMap single] Lifted $N_LIFTED / $N_INPUT variants ($N_DROPPED dropped)" >> {log} 2>&1
                    awk '{{print $4}}' {output.tempDir}/study_hg38.bed > {output.tempDir}/lifted_snps.txt
                    awk '{{print $4, $3}}' {output.tempDir}/study_hg38.bed > {output.tempDir}/new_pos.txt
                    awk '{{gsub(/^chr/,"",$1); print $4, $1}}' {output.tempDir}/study_hg38.bed > {output.tempDir}/new_chr.txt

                    plink2 --pfile {params.input_prefix} --extract {output.tempDir}/lifted_snps.txt --make-pgen --out {output.tempDir}/step1
                    plink2 --pfile {output.tempDir}/step1 --update-map {output.tempDir}/new_pos.txt --make-pgen --out {output.tempDir}/step2
                    plink2 --pfile {output.tempDir}/step2 --update-chr {output.tempDir}/new_chr.txt --sort-vars --set-all-var-ids 'chr@:#:$r:$a' --make-pgen --out {params.output_prefix}
                """)


rule crossmapLdPrunedToB38:
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/crossmap:latest"
    conda:
        "../../envs/crossmap.yml"
    envmodules: *([config.get("plink_module")] if config.get("plink_module") else [])
    threads: 2
    resources:
        nodes=1,
        mem_mb=8000,
        runtime=30,
    input:
        LDpgen=OUT_DIR / "{subset}" / "f1.ldpruned.pgen",
        LDpvar=OUT_DIR / "{subset}" / "f1.ldpruned.pvar",
        LDpsam=OUT_DIR / "{subset}" / "f1.ldpruned.psam",
    output:
        LDpgen=OUT_DIR / "{subset}" / "f1.b38.ldpruned.pgen",
        LDpvar=OUT_DIR / "{subset}" / "f1.b38.ldpruned.pvar",
        LDpsam=OUT_DIR / "{subset}" / "f1.b38.ldpruned.psam",
        tempDir=temp(directory(OUT_DIR / "{subset}" / "intermediates" / "crossmap_ld")),
    params:
        ld_prefix=lambda wildcards, input: input.LDpgen[:-5],
        output_prefix=lambda wildcards: str(OUT_DIR / wildcards.subset / "f1.b38.ldpruned"),
    run:
        _build = config.get("build", "GRCh38")
        if _build == "GRCh38":
            shell("""
                mkdir -p {output.tempDir}
                plink2 --pfile {params.ld_prefix} \
                       --make-pgen \
                       --out {params.output_prefix}
            """)
        else:
            shell("""
                mkdir -p {output.tempDir}
                awk '$1 !~ /^#/ {{print $3}}' {params.ld_prefix}.pvar > {output.tempDir}/ldpruned_vars.txt
                plink2 --pfile {params.ld_prefix} \
                       --extract {output.tempDir}/ldpruned_vars.txt \
                       --make-pgen \
                       --out {params.output_prefix}
            """)