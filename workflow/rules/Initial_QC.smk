rule initialFilter :
    container: "oras://ghcr.io/coffm049/gdcgnomicsqc/plink:latest"
    conda: "../../envs/ancNreport.yml"
    threads: 8
    resources:
        nodes = 1,
        mem_mb = 16000,
        runtime = 60,
    output:
        bed = OUT_DIR / "{subset}" / "initialFilter.bed",
        bim = OUT_DIR / "{subset}" / "initialFilter.bim",
        fam = OUT_DIR / "{subset}" / "initialFilter.fam",
        LDbed = OUT_DIR / "{subset}" / "initialFilter.LDpruned.bed",
        LDbim = OUT_DIR / "{subset}" / "initialFilter.LDpruned.bim",
        LDfam = OUT_DIR / "{subset}" / "initialFilter.LDpruned.fam",
        tempDir  = temp(directory(OUT_DIR / "{subset}" / "intermediates" / "initial_filter")),
        smiss = OUT_DIR / "{subset}" / "initial.smiss",
        vmiss = OUT_DIR / "{subset}" / "initial.vmiss",
        smissIMG = report(OUT_DIR / "{subset}" / "figures" / "smiss.svg", caption = "../../report/smiss.rst", category = "Quality Control"),
        vmissIMG = report(OUT_DIR / "{subset}" / "figures" / "vmiss.svg", caption = "../../report/vmiss.rst", category = "Quality Control"),
    input:
        pgen = expand(OUT_DIR / "{{subset}}" / "initialFilter_{CHR}.pgen", CHR = CHROMOSOMES), 
        psam = expand(OUT_DIR / "{{subset}}" / "initialFilter_{CHR}.psam", CHR = CHROMOSOMES), 
        pvar = expand(OUT_DIR / "{{subset}}" / "initialFilter_{CHR}.pvar", CHR = CHROMOSOMES), 
        ancestries = get_ancestry_file # Snakemake evaluates this per wildcard
    params:
        output_prefix = lambda wildcards, output: output.bed[:-4],
    shell: """
    
    # MERGE chromosomes
    # 1. Create/Clear the merge list file

    mkdir -p {output.tempDir}
    > {output.tempDir}/mergelist.txt

    # 2. Extract prefixes and write to the list
    for f in {input.pgen}; do
        echo "${{f%.pgen}}" >> {output.tempDir}/mergelist.txt
    done

    plink2 --pmerge-list {output.tempDir}/mergelist.txt \
           --make-bed \
           --missing \
           --out {output.tempDir}/intermediate_0


    if [ "{wildcards.subset}" == "full" ]; then
         echo "Processing full dataset without subsetting..."
         
         bash scripts/initialFilter.sh {output.tempDir}/intermediate_0 {params.output_prefix} {threads}  {output.tempDir}
                
    else
         echo "Subsetting data for ancestry: {wildcards.subset}..."
         
         bash scripts/initialFilter.sh {output.tempDir}/intermediate_0 {params.output_prefix} {threads} {output.tempDir}
    fi
    mv {output.tempDir}/intermediate_0.vmiss {output.vmiss}
    mv {output.tempDir}/intermediate_0.smiss {output.smiss}
    
    Rscript scripts/plotMissingness.R {output.smiss} {output.vmiss} {output.smissIMG} {output.vmissIMG} 
    """
