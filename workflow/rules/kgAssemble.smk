rule kgAssemble:
    log:
        OUT_DIR / "logs" / "kgAssemble.log",
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
    conda:
        "../../envs/ancNreport.yml"
    threads: 8
    resources:
        nodes=1,
        mem_mb=32000,
        runtime=30,
    output:
        tempDir=temp(directory(REF / "intermediates")),
        highcovPgen=protected(
            REF / "1000G_highcoverage" / "1000G_highCoveragephased.pgen"
        ),
        highcovPvar=protected(
            REF / "1000G_highcoverage" / "1000G_highCoveragephased.pvar"
        ),
        highcovPsam=protected(
            REF / "1000G_highcoverage" / "1000G_highCoveragephased.psam"
        ),
        ldPgen=protected(
            REF / "1000G_highcoverage" / "1000G_highCoveragephased.pruned.pgen"
        ),
        ldPvar=protected(
            REF / "1000G_highcoverage" / "1000G_highCoveragephased.pruned.pvar"
        ),
        ldPsam=protected(
            REF / "1000G_highcoverage" / "1000G_highCoveragephased.pruned.psam"
        ),
        fastafai=protected(REF / "Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai"),
    input:
        vcf=expand(
            REF
            / "1000G_highcoverage"
            / "1kGP_high_coverage_Illumina.chr{CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz",
            CHR=CHROMOSOMES,
        ),
        fasta=REF / "Homo_sapiens.GRCh38.dna.primary_assembly.fa",
    params:
        highcovPgen=REF / "1000G_highcoverage" / "1000G_highCoveragephased",
        ldPgen=REF / "1000G_highcoverage" / "1000G_highCoveragephased.pruned",
    shell:
        """
        samtools faidx {input.fasta}
        mkdir -p {output.tempDir}
        rm -f {output.tempDir}/mergelist.txt
        for f in {input.vcf}; do
            FILE_NAME=$(basename ${{f}})

            plink2 --vcf ${{f}} \
                --make-pgen \
                --maf 0.05 \
                --threads {threads} \
                --snps-only 'just-acgt' \
                --rm-dup force-first \
                --out {output.tempDir}/${{FILE_NAME%.vcf.gz}}
            echo "{output.tempDir}/${{FILE_NAME%.vcf.gz}}" >> {output.tempDir}/mergelist.txt
        done

        plink2 --pmerge-list {output.tempDir}/mergelist.txt \
            --make-pgen \
            --threads {threads} \
            --out {output.tempDir}/full
        

        plink2 --pfile {output.tempDir}/full \
            --make-pgen \
            --threads {threads} \
            --ref-from-fa --fa {input.fasta} \
            --out {output.tempDir}/full2
        plink2 --pfile {output.tempDir}/full2 \
            --make-pgen \
            --threads {threads} \
            --set-all-var-ids 'chr@:#:$r:$a' \
            --out {output.tempDir}/full3

        plink2 --pfile {output.tempDir}/full3 \
               --indep-pairphase 1000kb 1 0.1 \
               --threads {threads} \
               --out {output.tempDir}/full3
        
        plink2 --pfile {output.tempDir}/full3 \
               --threads {threads} \
               --extract {output.tempDir}/full3.prune.in \
               --king-cutoff 0.0884 \
               --out {output.tempDir}/full3


        # only unrelated in reference sets
        plink2 --pfile {output.tempDir}/full3 \
            --threads {threads} \
            --keep {output.tempDir}/full3.king.cutoff.in.id \
            --make-pgen \
            --out {params.highcovPgen}


        plink2 --pfile {output.tempDir}/full3 \
            --threads {threads} \
            --make-pgen \
            --set-all-var-ids 'chr@:#:$r:$a' \
            --keep {output.tempDir}/full3.king.cutoff.in.id \
            --extract {output.tempDir}/full3.prune.in \
            --out {params.ldPgen}
        """
