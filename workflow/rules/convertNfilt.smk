rule convertNfilt :
    container: "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
    conda: "../../envs/ancNreport.yml"
    threads: 4
    resources:
        nodes = 1,
        mem_mb = 32000,
        runtime = 60,
    output:
        pgen = OUT_DIR / "{subset}" / "initialFilter_{CHR}.pgen",
        pvar = OUT_DIR / "{subset}" / "initialFilter_{CHR}.pvar",
        psam = OUT_DIR / "{subset}" / "initialFilter_{CHR}.psam",
        tempDir  = temp(directory(OUT_DIR / "{subset}" / "{CHR}" / "intermediates" / "initial_filter")),
        smiss = OUT_DIR / "{subset}" / "initial_{CHR}.smiss",
        vmiss = OUT_DIR / "{subset}" / "initial_{CHR}.vmiss"
    input:
        vcf = config['vcf_template'],
        ancestries = get_ancestry_file, # Snakemake evaluates this per wildcard
        crossmap = REF / "CrossMap" / "hg19ToHg38.over.chain.gz",
        gr38fasta = REF / "Homo_sapiens.GRCh38.dna.primary_assembly.fa",
    params:
        thin = config['thin'],
        input_prefix = lambda wildcards, input: input.vcf[:-7],
        output_prefix = lambda wildcards, output: output.pgen[:-5],
        liftover = True
    shell: """
    
    mkdir -p {output.tempDir}
    # if [[ "{params.liftover}" == "True" ]]; then
    #     for i in {{1..22}} X Y MT; do echo "$i chr$i"; done > {output.tempDir}/chr_map.txt
    #     
    #     # Use it in the command
    #     bcftools annotate --rename-chrs {output.tempDir}/chr_map.txt {input.vcf} -Ou \
    #     | bcftools +liftover -- \
    #       -c {input.crossmap} \
    #       -f {input.gr38fasta} \
    #       -o {output.tempDir}/intermediate_00.vcf.gz \
    #       --reject {output.tempDir}/rejected_variants.vcf.gz
    # else 
    #     ln -sf $(realpath {input.vcf}) {output.tempDir}/intermediate_00.vcf.gz
    # fi

    if [[ "{wildcards.subset}" == "full" && "{params.thin}" == "True" ]]; then
        echo "Processing full dataset without subsetting..."
        plink2 --vcf {input.vcf} \
               --make-pgen \
               --rm-dup force-first \
               --missing \
               --thin-indiv-count 5000 \
               --thin-count 20000 \
               --threads {threads} \
               --seed 1 \
               --memory {resources.mem_mb} \
               --out {output.tempDir}/intermediate_0
    elif [[ "{wildcards.subset}" == "full" && "{params.thin}" != "True" ]]; then
        plink2 --vcf {input.vcf} \
               --make-pgen \
               --rm-dup force-first \
               --threads {threads} \
               --missing \
               --memory {resources.mem_mb} \
               --out {output.tempDir}/intermediate_0
    elif [[ "{wildcards.subset}" != "full" && "{params.thin}" == "True" ]]; then
        plink2 --vcf {input.vcf} \
               --make-pgen \
               --rm-dup force-first \
               --threads {threads} \
               --missing \
               --covar {input.ancestries} \
               --keep-if 'pc_label == {wildcards.subset}' \
               --thin-indiv-count 5000 \
               --thin-count 20000 \
               --seed 1 \
               --memory {resources.mem_mb} \
               --out {output.tempDir}/intermediate_0
    elif [[ "{wildcards.subset}" != "full" && "{params.thin}" != "True" ]]; then
        plink2 --vcf {input.vcf} \
               --make-pgen \
               --missing \
               --threads {threads} \
               --rm-dup force-first \
               --memory {resources.mem_mb} \
               --covar {input.ancestries} \
               --keep-if 'pc_label == {wildcards.subset}' \
               --out {output.tempDir}/intermediate_0
    fi

    plink2 --pfile {output.tempDir}/intermediate_0 \
           --make-pgen \
           --geno 0.1 \
           --threads {threads} \
           --snps-only 'just-acgt' \
           --output-chr 26 \
           --sort-vars \
           --out {params.output_prefix}


    mv {output.tempDir}/intermediate_0.vmiss {output.vmiss}
    mv {output.tempDir}/intermediate_0.smiss {output.smiss}
    """
