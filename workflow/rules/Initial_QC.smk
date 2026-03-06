def get_ancestry_file(wildcards):
    """
    If the subset is 'full', return an empty list (no dependency).
    Otherwise, return the path to the estimated ancestry file.
    """
    if wildcards.subset == "full":
        return []
    else:
        # This forces Snakemake to wait for the estimation rule to finish
        # OUT_DIR is globally defined
        return OUT_DIR / "02-globalAncestry" / "latentDistantRelatedness.csv",


rule initialFilter :
    container: "oras://ghcr.io/coffm049/gdcgnomicsqc/plink:latest"
    conda: "../../envs/ancNreport.yml"
    threads: 8
    resources:
        nodes = 1,
        mem_mb = 16000,
        runtime = 60,
    input:
        vcf = f"{config['vcf_template']}.vcf.gz",
        ancestries = get_ancestry_file # Snakemake evaluates this per wildcard
    output:
        bed = OUT_DIR / "{subset}" / "initialFilter_{CHR}.bed",
        bim = OUT_DIR / "{subset}" / "initialFilter_{CHR}.bim",
        fam = OUT_DIR / "{subset}" / "initialFilter_{CHR}.fam",
        LDbed = OUT_DIR / "{subset}" / "initialFilter_{CHR}.LDpruned.bed",
        LDbim = OUT_DIR / "{subset}" / "initialFilter_{CHR}.LDpruned.bim",
        LDfam = OUT_DIR / "{subset}" / "initialFilter_{CHR}.LDpruned.fam",
        tempDir  = temp(directory(OUT_DIR / "{subset}" / "{CHR}" / "intermediates" / "initial_filter")),
        smiss = OUT_DIR / "{subset}" / "initial_{CHR}.smiss",
        vmiss = OUT_DIR / "{subset}" / "initial_{CHR}.vmiss"
    params:
        thin = config['thin'],
        input_prefix = lambda wildcards, input: input.vcf[:-7],
        output_prefix = lambda wildcards, output: output.bed[:-4],
    shell: """
    mkdir -p {output.tempDir}

    if [ "{wildcards.subset}" == "full" ]; then
         echo "Processing full dataset without subsetting..."
         
         plink2 --vcf {input.vcf} \
                --make-bed \
                --missing \
                --out {output.tempDir}/intermediate_0
         bash scripts/initialFilter.sh {output.tempDir}/intermediate_0 {params.output_prefix} {threads} {params.thin} {output.tempDir}
                
    else
         echo "Subsetting data for ancestry: {wildcards.subset}..."
         
         # The quotes around "Ancestry == {wildcards.subset}" ensure Bash 
         # passes the entire condition safely to PLINK2.
         # (Note: Ensure 'Ancestry' matches your actual column header)
         plink2 --vcf {input.vcf} \
                --pheno {input.ancestries} \
                --keep-if 'pc_label == "'{wildcards.subset}'"'
                --make-bed \
                --missing \
                --out {output.tempDir}/intermediate_0
         bash scripts/initialFilter.sh {output.tempDir}/intermediate_0 {params.output_prefix} {threads} {params.thin} {output.tempDir}
    fi
    mv {output.tempDir}/intermediate_0.vmiss {output.vmiss}
    mv {output.tempDir}/intermediate_0.smiss {output.smiss}
    """
