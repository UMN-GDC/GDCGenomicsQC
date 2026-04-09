rule convertVcfToPlinkPerChromosome:
    log:
        OUT_DIR / "logs" / "convertVcfToPlinkPerChromosome_{subset}_{CHR}.log",
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
    conda:
        "../../envs/ancNreport.yml"
    threads: 4
    resources:
        nodes=1,
        mem_mb=32000,
        runtime=240,
    output:
        pgen=OUT_DIR / "{subset}" / "initialFilter_{CHR}.pgen",
        pvar=OUT_DIR / "{subset}" / "initialFilter_{CHR}.pvar",
        psam=OUT_DIR / "{subset}" / "initialFilter_{CHR}.psam",
        tempDir=temp(
            directory(
                OUT_DIR / "{subset}" / "{CHR}" / "intermediates" / "initial_filter"
            )
        ),
        smiss=OUT_DIR / "{subset}" / "initial_{CHR}.smiss",
        vmiss=OUT_DIR / "{subset}" / "initial_{CHR}.vmiss",
    input:
        vcf=config.get("vcf_template", ""),
        keep=get_ancestry_file,
        crossmap=REF / "CrossMap" / "hg19ToHg38.over.chain.gz",
        gr38fasta=REF / "Homo_sapiens.GRCh38.dna.primary_assembly.fa",
    params:
        thin=config.get("thin", False),
        # Parameters for plink2 filtering
        min_mach_r2=config.get("convertNfilt", {}).get("info_r2_min"),
        max_mach_r2=config.get("convertNfilt", {}).get("info_r2_max"),
        qual_min=config.get("convertNfilt", {}).get("qual_min"),
        output_prefix=lambda wildcards, output: output.pgen.replace(".pgen", ""),
    shell:
        """
    mkdir -p {output.tempDir}

    # Prepare plink2 filtering flags
    PLINK2_FILTERS=""
    if [ -n "{params.min_mach_r2}" ] && [ "{params.min_mach_r2}" != "None" ]; then
        if [ -n "{params.max_mach_r2}" ] && [ "{params.max_mach_r2}" != "None" ]; then
            PLINK2_FILTERS="$PLINK2_FILTERS --mach-r2-filter {params.min_mach_r2} {params.max_mach_r2}"
        fi
    fi
    if [ -n "{params.qual_min}" ] && [ "{params.qual_min}" != "None" ] && [ "{params.qual_min}" != "0" ]; then
        PLINK2_FILTERS="$PLINK2_FILTERS --var-min-qual {params.qual_min}"
    fi

    echo "Applying filters: $PLINK2_FILTERS"

    # Base PLINK2 command for conversion and filtering
    # Start constructing command
    CMD="plink2 --vcf {input.vcf} --make-pgen --rm-dup force-first --missing --threads {threads} --memory {resources.mem_mb} --out {output.tempDir}/intermediate_0 $PLINK2_FILTERS"

    if [[ "{wildcards.subset}" != "full" ]]; then
        CMD="$CMD --keep {input.keep}"
    fi

    if [[ "{params.thin}" == "True" ]]; then
        if [[ "{wildcards.subset}" == "full" ]]; then
            CMD="$CMD --thin-indiv 0.1 --thin-count 100000 --seed 1"
        else
            CMD="$CMD --thin-indiv-count 10000 --thin-count 100000 --seed 1"
        fi
    fi

    # Execute main conversion/filter command
    $CMD

    # Secondary processing
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
