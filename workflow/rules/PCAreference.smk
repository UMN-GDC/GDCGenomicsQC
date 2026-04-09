rule PCAreference:
    log:
        OUT_DIR / "logs" / "PCAreference.log",
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
    conda:
        "../../envs/ancNreport.yml"
    threads: 8
    resources:
        nodes=1,
        mem_mb=32000,
        runtime=2880,
    input:
        pgen=OUT_DIR / "full" / "initialFilter.pgen",
        pvar=OUT_DIR / "full" / "initialFilter.pvar",
        psam=OUT_DIR / "full" / "initialFilter.psam",
        ldPgen=REF / "1000G_highcoverage" / "1000G_highCoveragephased.pruned.pgen",
        ldPvar=REF / "1000G_highcoverage" / "1000G_highCoveragephased.pruned.pvar",
        ldPsam=REF / "1000G_highcoverage" / "1000G_highCoveragephased.pruned.psam",
    output:
        # List all files that PLINK will actually create
        eigen=OUT_DIR / "01-globalAncestry" / "ref.eigenvec",
        projected=OUT_DIR / "01-globalAncestry" / "sampleRefPCscores.sscore",
        projectedref=OUT_DIR / "01-globalAncestry" / "refRefPCscores.sscore",
        tempDir=temp(directory(OUT_DIR / "01-globalAncestry" / "intermediates")),
    params:
        method=config.get("relatedness", {}).get("method", "king"),
        grm=config.get("relatedness", {}).get("method", "king"),
        input_prefix=OUT_DIR / "full" / "initialFilter",
        dir=str(OUT_DIR / "01-globalAncestry"),
        ref=REF / "1000G_highcoverage" / "1000G_highCoveragephased.pruned",
    shell:
        """
        echo "PCA: "
        mkdir -p {output.tempDir}
        
        # Filter pruned ref panel to shared SNPs and LDprune
        plink2 --pfile {params.input_prefix} --write-snplist \
            --maf 0.05 \
            --threads {threads} \
            --out {params.dir}/intermediates/study_snps
        

        # calculate ld on ref intersection with sample with common frequency
        # compute PCA refrence only
        plink2 --pfile {params.ref} \
               --freq counts \
               --threads {threads} \
               --extract {params.dir}/intermediates/study_snps.snplist \
               --pca approx allele-wts vcols=chrom,ref,alt \
               --out {params.dir}/ref \
               --allow-no-sex
        
        echo "Project sample onto the reference PCs."
        plink2 --pfile {params.input_prefix} \
               --read-freq {params.dir}/ref.acount \
               --score {params.dir}/ref.eigenvec.allele 2 5 header-read \
               --score-col-nums 6-15 \
               --out {params.dir}/sampleRefPCscores
        echo "Project ref onto the reference PCs."
        
        plink2 --pfile {params.ref} \
               --read-freq {params.dir}/ref.acount \
               --score {params.dir}/ref.eigenvec.allele 2 5 header-read \
               --score-col-nums 6-15 \
               --out {params.dir}/refRefPCscores
        """
