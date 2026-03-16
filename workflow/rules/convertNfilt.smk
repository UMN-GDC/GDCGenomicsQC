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
        ancestries = get_ancestry_file # Snakemake evaluates this per wildcard
    params:
        thin = config['thin'],
        input_prefix = lambda wildcards, input: input.vcf[:-7],
        output_prefix = lambda wildcards, output: output.pgen[:-5],
    shell: """
    mkdir -p {output.tempDir}

    if [[ "{wildcards.subset}" == "full" && "{params.thin}" == "True" ]]; then
        echo "Processing full dataset without subsetting..."
        plink2 --vcf {input.vcf} \
               --make-pgen \
               --missing \
               --thin-indiv-count 5000 \
               --thin-count 5000 \
               --seed 1 \
               --memory {resources.mem_mb} \
               --out {output.tempDir}/intermediate_0
    elif [[ "{wildcards.subset}" == "full" && "{params.thin}" != "True" ]]; then
        plink2 --vcf {input.vcf} \
               --make-pgen \
               --missing \
               --memory {resources.mem_mb} \
               --out {output.tempDir}/intermediate_0
    elif [[ "{wildcards.subset}" != "full" && "{params.thin}" == "True" ]]; then
        plink2 --vcf {input.vcf} \
               --make-pgen \
               --missing \
               --covar {input.ancestries} \
               --keep-if 'pc_label == {wildcards.subset}' \
               --thin-indiv-count 5000 \
               --thin-count 5000 \
               --seed 1 \
               --memory {resources.mem_mb} \
               --out {output.tempDir}/intermediate_0
    elif [[ "{wildcards.subset}" != "full" && "{params.thin}" != "True" ]]; then
        plink2 --vcf {input.vcf} \
               --make-pgen \
               --missing \
               --memory {resources.mem_mb} \
               --covar {input.ancestries} \
               --keep-if 'pc_label == {wildcards.subset}' \
               --out {output.tempDir}/intermediate_0
    fi

    plink2 --pfile {output.tempDir}/intermediate_0 \
           --make-pgen \
           --geno 0.1 \
           --memory {resources.mem_mb} \
           --out {params.output_prefix}

    mv {output.tempDir}/intermediate_0.vmiss {output.vmiss}
    mv {output.tempDir}/intermediate_0.smiss {output.smiss}
    """
