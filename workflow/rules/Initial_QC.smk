rule mergeChromosomesAndFilter:
    log:
        OUT_DIR / "logs" / "mergeChromosomesAndFilter_{subset}.log",
    container:
        "docker://gfanz/plink2:latest"
    conda:
        "../../envs/ancNreport.yml"
    envmodules: use("plink_module")
    threads: 8
    resources:
        nodes=1,
        mem_mb=128000,
        runtime=240,
    output:
        pgen=OUT_DIR / "{subset}" / "f1.pgen",
        pvar=OUT_DIR / "{subset}" / "f1.pvar",
        psam=OUT_DIR / "{subset}" / "f1.psam",
        LDpgen=OUT_DIR / "{subset}" / "f1.ldpruned.pgen",
        LDpvar=OUT_DIR / "{subset}" / "f1.ldpruned.pvar",
        LDpsam=OUT_DIR / "{subset}" / "f1.ldpruned.psam",
        tempDir=temp(
            directory(OUT_DIR / "{subset}" / "intermediates" / "initial_filter")
        ),
        smiss=OUT_DIR / "{subset}" / "initial.smiss",
        vmiss=OUT_DIR / "{subset}" / "initial.vmiss",
    input:
        fasta=ancient(REF / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"),
        pgen=expand(
            OUT_DIR / "{{subset}}" / "f1_{CHR}.pgen", CHR=CHROMOSOMES
        ),
        pvar=expand(
            OUT_DIR / "{{subset}}" / "f1_{CHR}.pvar", CHR=CHROMOSOMES
        ),
        psam=expand(
            OUT_DIR / "{{subset}}" / "f1_{CHR}.psam", CHR=CHROMOSOMES
        ),
        keep=get_ancestry_file,
    params:
        output_prefix=lambda wildcards, output: output.pgen[:-5],
        scripts_dir=SCRIPTS_DIR,
    shell:
        """
    
    # MERGE chromosomes
    # 1. Create/Clear the merge list file

    mkdir -p {output.tempDir}
    > {output.tempDir}/mergelist.txt

    # 2. Extract prefixes and write to the list
    for f in {input.pgen}; do
        echo "${{f%.pgen}}" >> {output.tempDir}/mergelist.txt
    done

    plink2 --pmerge-list {output.tempDir}/mergelist.txt \
           --threads {threads} \
           --make-pgen \
           --missing \
           --out {output.tempDir}/intermediate_0
    plink2 --pfile {output.tempDir}/intermediate_0 \
           --threads {threads} \
           --fa {input.fasta} --ref-from-fa force \
           --out {output.tempDir}/intermediate_1 \
           --make-pgen
    plink2 --pfile {output.tempDir}/intermediate_1 \
           --threads {threads} \
           --set-all-var-ids 'chr@:#:$r:$a' \
           --out {output.tempDir}/intermediate_2 \
           --make-pgen


    if [ "{wildcards.subset}" == "full" ]; then
         echo "Processing full dataset without subsetting..."
         
         bash {params.scripts_dir}/initialFilter.sh {output.tempDir}/intermediate_2 {params.output_prefix} {threads}  {output.tempDir}
                 
    else
         echo "Subsetting data for ancestry: {wildcards.subset}..."
         
         bash {params.scripts_dir}/initialFilter.sh {output.tempDir}/intermediate_2 {params.output_prefix} {threads} {output.tempDir}
    fi
    mv {output.tempDir}/intermediate_0.vmiss {output.vmiss}
    mv {output.tempDir}/intermediate_0.smiss {output.smiss}
    for ext in pgen pvar psam; do
        mv {params.output_prefix}.LDpruned.$ext {params.output_prefix}.ldpruned.$ext
    done
    """
