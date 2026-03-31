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
        keep = get_ancestry_file,
        crossmap = REF / "CrossMap" / "hg19ToHg38.over.chain.gz",
        gr38fasta = REF / "Homo_sapiens.GRCh38.dna.primary_assembly.fa",
    params:
        thin = config['thin'],
        input_prefix = lambda wildcards, input: input.vcf[:-7],
        output_prefix = lambda wildcards, output: output.pgen[:-5],
        liftoover = True,
        info_r2_min = config.get('convertNfilt', {}).get('info_r2_min'),
        filter_pass = config.get('convertNfilt', {}).get('filter_pass', True),
        qual_min = config.get('convertNfilt', {}).get('qual_min'),
    shell: """
    
    mkdir -p {output.tempDir}
    
    BCFTOOLS_FILTER=""
    {params.filter_pass:+BCFTOOLS_FILTER="$BCFTOOLS_FILTER -i 'FILTER==\"PASS\"'"}
    {params.qual_min:+BCFTOOLS_FILTER="$BCFTOOLS_FILTER -i 'QUAL>={params.qual_min}'"}
    {params.info_r2_min:+BCFTOOLS_FILTER="$BCFTOOLS_FILTER -i 'INFO/R2>={params.info_r2_min}'"}
    
    # if [[ "{params.liftoover}" == "True" ]]; then
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

    if [ -n "$BCFTOOLS_FILTER" ]; then
        echo "Pre-filtering VCF with bcftools: $BCFTOOLS_FILTER"
        bcftools view {input.vcf} $BCFTOOLS_FILTER -Oz -o {output.tempDir}/filtered.vcf.gz
        VCF_INPUT={output.tempDir}/filtered.vcf.gz
    else
        VCF_INPUT={input.vcf}
    fi

    if [[ "{wildcards.subset}" == "full" && "{params.thin}" == "True" ]]; then
        echo "Processing full dataset without subsetting..."
        plink2 --vcf $VCF_INPUT \
               --make-pgen \
               --rm-dup force-first \
               --missing \
               --thin-indiv 0.1 \
               --thin-count 100000 \
               --threads {threads} \
               --seed 1 \
               --memory {resources.mem_mb} \
               --out {output.tempDir}/intermediate_0
    elif [[ "{wildcards.subset}" == "full" && "{params.thin}" != "True" ]]; then
        plink2 --vcf $VCF_INPUT \
               --make-pgen \
               --rm-dup force-first \
               --threads {threads} \
               --missing \
               --memory {resources.mem_mb} \
               --out {output.tempDir}/intermediate_0
    elif [[ "{wildcards.subset}" != "full" && "{params.thin}" == "True" ]]; then
        plink2 --vcf $VCF_INPUT \
               --make-pgen \
               --rm-dup force-first \
               --threads {threads} \
               --missing \
               --keep {input.keep} \
               --thin-indiv-count 10000 \
               --thin-count 100000 \
               --seed 1 \
               --memory {resources.mem_mb} \
               --out {output.tempDir}/intermediate_0
    elif [[ "{wildcards.subset}" != "full" && "{params.thin}" != "True" ]]; then
        plink2 --vcf $VCF_INPUT \
               --make-pgen \
               --missing \
               --threads {threads} \
               --rm-dup force-first \
               --memory {resources.mem_mb} \
               --keep {input.keep} \
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
