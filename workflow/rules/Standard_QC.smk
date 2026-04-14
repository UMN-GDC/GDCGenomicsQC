rule applyStandardQualityControl:
    log:
        OUT_DIR / "logs" / "applyStandardQualityControl_{subset}.log",
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
    conda:
        "../../envs/ancNreport.yml"
    threads: 8
    resources:
        nodes=1,
        mem_mb=32000,
        runtime=60,
    input:
        pgen=OUT_DIR / "{subset}" / "initialFilter.pgen",
        pvar=OUT_DIR / "{subset}" / "initialFilter.pvar",
        psam=OUT_DIR / "{subset}" / "initialFilter.psam",
        LDpgen=OUT_DIR / "{subset}" / "initialFilter.LDpruned.pgen",
        LDpvar=OUT_DIR / "{subset}" / "initialFilter.LDpruned.pvar",
        LDpsam=OUT_DIR / "{subset}" / "initialFilter.LDpruned.psam",
    output:
        pgen=OUT_DIR / "{subset}" / "standardFilter.pgen",
        pvar=OUT_DIR / "{subset}" / "standardFilter.pvar",
        psam=OUT_DIR / "{subset}" / "standardFilter.psam",
        LDpgen=OUT_DIR / "{subset}" / "standardFilter.LDpruned.pgen",
        LDpvar=OUT_DIR / "{subset}" / "standardFilter.LDpruned.pvar",
        LDpsam=OUT_DIR / "{subset}" / "standardFilter.LDpruned.psam",
        tempDir=temp(
            directory(OUT_DIR / "{subset}" / "intermediates" / "standard_filter")
        ),
    params:
        ref=config.get("REF", "/path/to/ref"),
        output_dir=lambda wildcards, input: OUT_DIR / wildcards.subset,
        sex_check=config.get("SEX_CHECK", False),
        input_prefix=lambda wildcards, input: input.LDpgen[:-5],
        relatedness=config.get("relatedness", {}).get("method", "king"),
    shell:
        """
        echo "Standard QC: Variants and samples filtering"
        echo "Data subset: {wildcards.subset}"
        mkdir -p {output.tempDir}

        if [[ "{params.sex_check}" == "True" ]]; then
          echo "Performing Sex check"
          plink2 --pfile {params.input_prefix} --check-sex --out {params.input_prefix} --threads {threads}
          grep 'PROBLEM' {params.input_prefix}.sexcheck | awk '{{print $1,$2}}' > {params.output_dir}/sex_discrepancy.txt
          plink2 --pfile {params.input_prefix} --remove {params.output_dir}/sex_discrepancy.txt --make-pgen --out {output.tempDir}/pastSex --threads {threads}
        else
          echo "Ignoring Sex check"
          mv {input.LDpgen} {output.tempDir}/pastSex.pgen
          mv {input.LDpvar} {output.tempDir}/pastSex.pvar
          mv {input.LDpsam} {output.tempDir}/pastSex.psam
        fi
        bash scripts/filterStandard.sh {output.tempDir}/pastSex {params.output_dir} {threads}
        """
