rule PCAreference:
    container: "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
    conda: "../../envs/ancNreport.yml"
    threads: 8
    resources:
        nodes = 1,
        mem_mb = 32000,
        runtime = 2880,
    input:
        bed = OUT_DIR / "full" / "initialFilter.bed",
        bim = OUT_DIR / "full" / "initialFilter.bim",
        fam = OUT_DIR / "full" / "initialFilter.fam",
    output:
        # List all files that PLINK will actually create
        eigen = OUT_DIR / "01-globalAncestry" / "ref.eigenvec",
        projected = OUT_DIR / "01-globalAncestry" / "sampleRefPCscores.sscore",
        projectedref = OUT_DIR / "01-globalAncestry" / "refRefPCscores.sscore",
        tempDir  = temp(directory(OUT_DIR / "01-globalAncestry" / "intermediates")),
    params:
        method = config['relatedness']["method"],
        grm = config['relatedness']["method"],
        input_prefix = OUT_DIR / "full" / "initialFilter",
        dir = directory(OUT_DIR / "01-globalAncestry"),
        ref= REF / "1000G_GRCh38" / "1000G.ensembl.105.with.rsid.gender"
    shell: 
        """
        echo "PCA: "
        mkdir -p {output.tempDir}
        
        # Filter ref panel to shared SNPs and LDprune
        plink2 --bfile {params.input_prefix} --write-snplist \
            --maf 0.05 \
            --rm-dup force-first --out {params.dir}/intermediates/study_snps
        plink2 --bfile {params.ref} \
               --extract {params.dir}/intermediates/study_snps.snplist \
               --maf 0.05 \
               --indep-pairwise 200 50 0.1 \
               --out {output.tempDir}/ref_filtered \
               --make-bed
        
        # compute PCA refrence only
        if [[ "{params.grm}" == "True" ]] ; then
          plink2 --bfile {output.tempDir}/ref_filtered \
                 --freq counts \
                 --extract {output.tempDir}/ref_filtered.prune.in \
                 --pca approx allele-wts vcols=chrom,ref,alt \
                 --out {params.dir}/ref \
                 --make-grm-bin \
                 --allow-no-sex
        else 
          plink2 --bfile {output.tempDir}/ref_filtered \
                 --freq counts \
                 --extract {output.tempDir}/ref_filtered.prune.in \
                 --pca approx allele-wts vcols=chrom,ref,alt \
                 --out {params.dir}/ref \
                 --allow-no-sex
        fi
        
        echo "Project sample onto the reference PCs."
        plink2 --bfile {params.input_prefix} \
               --read-freq {params.dir}/ref.acount \
               --score {params.dir}/ref.eigenvec.allele 2 5 header-read \
               --score-col-nums 6-15 \
               --out {params.dir}/sampleRefPCscores
        echo "Project ref onto the reference PCs."
        
        plink2 --bfile {output.tempDir}/ref_filtered \
               --read-freq {params.dir}/ref.acount \
               --score {params.dir}/ref.eigenvec.allele 2 5 header-read \
               --score-col-nums 6-15 \
               --out {params.dir}/refRefPCscores
        """
