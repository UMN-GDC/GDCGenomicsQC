rule Standard_QC:
    container: "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
    conda: "../../envs/ancNreport.yml"
    threads: 8
    resources:
        # nodes=1 is usually the default, but can be specified if needed
        nodes = 1,
        # mem=32GB translated to MB
        mem_mb = 32000,
        runtime =60,
    input:
        bed = OUT_DIR / "{subset}" / "initialFilter.bed",
        bim = OUT_DIR / "{subset}" / "initialFilter.bim",
        fam = OUT_DIR / "{subset}" / "initialFilter.fam",
        LDbed = OUT_DIR / "{subset}" / "initialFilter.LDpruned.bed",
        LDbim = OUT_DIR / "{subset}" / "initialFilter.LDpruned.bim",
        LDfam = OUT_DIR / "{subset}" / "initialFilter.LDpruned.fam",
    output:
        bed =   OUT_DIR / "{subset}" / "standardFilter.bed",
        bim =   OUT_DIR / "{subset}" / "standardFilter.bim",
        fam =   OUT_DIR / "{subset}" / "standardFilter.fam",
        LDbed = OUT_DIR / "{subset}" / "standardFilter.LDpruned.bed",
        LDbim = OUT_DIR / "{subset}" / "standardFilter.LDpruned.bim",
        LDfam = OUT_DIR / "{subset}" / "standardFilter.LDpruned.fam",
        tempDir  = temp(directory(OUT_DIR / "{subset}" / "intermediates" / "standard_filter"))
    params:
        ref= config["REF"],
        output_dir = lambda wildcards, input: OUT_DIR / wildcards.subset,
        sex_check = config['SEX_CHECK'],
        input_prefix = lambda wildcards, input: input.LDbed[:-4],
        relatedness = config['relatedness']['method']
    shell: """
        echo "Standard QC: Variants and samples filtering"
        echo "Data subset: {wildcards.subset}"
        mkdir -p {output.tempDir}

        if [[ "{params.sex_check}" == "True" ]]; then
          echo "Performing Sex check"
          plink --bfile {params.input_prefix} --check-sex --out {params.input_prefix}
          grep 'PROBLEM' {params.input_prefix}.sexcheck | awk '{{print $1,$2}}' > {params.output_dir}/sex_discrepancy.txt
          plink --bfile {params.input_prefix} --remove {params.output_dir}/sex_discrepancy.txt --make-bed --out {output.tempDir}/pastSex
        else
          echo "Ignoring Sex check"
          mv {input.LDbed} {output.tempDir}/pastSex.bed
          mv {input.LDbim} {output.tempDir}/pastSex.bim
          mv {input.LDfam} {output.tempDir}/pastSex.fam
        fi
        bash scripts/filterStandard.sh {output.tempDir}/pastSex {params.output_dir} {threads}
        """
