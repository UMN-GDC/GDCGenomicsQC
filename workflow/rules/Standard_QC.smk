if INPUT_IS_PER_CHROMOSOME:
    rule applyStandardQualityControl:
        log:
            OUT_DIR / "logs" / "applyStandardQualityControl_{subset}_{CHR}.log",
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
            pgen=OUT_DIR / "{subset}" / "initialFilter_{CHR}.pgen",
            pvar=OUT_DIR / "{subset}" / "initialFilter_{CHR}.pvar",
            psam=OUT_DIR / "{subset}" / "initialFilter_{CHR}.psam",
        output:
            pgen=OUT_DIR / "{subset}" / "standardFilter_{CHR}.pgen",
            pvar=OUT_DIR / "{subset}" / "standardFilter_{CHR}.pvar",
            psam=OUT_DIR / "{subset}" / "standardFilter_{CHR}.psam",
            tempDir=temp(
                directory(OUT_DIR / "{subset}" / "intermediates" / "standard_filter_{CHR}")
            ),
        params:
            ref=config.get("REF", "/path/to/ref"),
            output_dir=lambda wildcards, input: OUT_DIR / wildcards.subset,
            sex_check=config.get("SEX_CHECK", False),
            input_prefix=lambda wildcards, input: input.pgen[:-5],
            relatedness=config.get("relatedness", {}).get("method", "king"),
            scripts_dir=SCRIPTS_DIR,
        shell:
            """
            echo "Standard QC: Variants and samples filtering"
            echo "Data subset: {wildcards.subset}"
            echo "Chromosome: {wildcards.CHR}"
            mkdir -p {output.tempDir}

            if [[ "{params.sex_check}" == "True" ]]; then
              echo "Performing Sex check"
              plink2 --pfile {params.input_prefix} --check-sex --out {params.input_prefix} --threads {threads}
              grep 'PROBLEM' {params.input_prefix}.sexcheck | awk '{{print $1,$2}}' > {params.output_dir}/sex_discrepancy_{wildcards.CHR}.txt
              plink2 --pfile {params.input_prefix} --remove {params.output_dir}/sex_discrepancy_{wildcards.CHR}.txt --make-pgen --out {output.tempDir}/step1 --threads {threads}
            else
              echo "Ignoring Sex check"
              cp {input.pgen} {output.tempDir}/step1.pgen
              cp {input.pvar} {output.tempDir}/step1.pvar
              cp {input.psam} {output.tempDir}/step1.psam
            fi

            plink2 --pfile {output.tempDir}/step1 --maf 0.01 --make-pgen --out {output.tempDir}/step2 --threads {threads}

            plink2 --pfile {output.tempDir}/step2 --hwe 1e-6 --make-pgen --out {output.tempDir}/step3a --threads {threads}
            plink2 --pfile {output.tempDir}/step3a --hwe 1e-10 --make-pgen --out {output.tempDir}/step3 --threads {threads}

            plink2 --pfile {output.tempDir}/step3 --exclude scripts/inversion.txt --range --indep-pairwise 50 5 0.2 --out {output.tempDir}/indepSNP --threads {threads}
            plink2 --pfile {output.tempDir}/step3 --extract {output.tempDir}/indepSNP.prune.in --het --out {output.tempDir}/hetcheck --threads {threads}

            Rscript --no-save scripts/heterozygosity_outliers_list.R
            if [ -f {params.output_dir}/het_fail_ind.txt ]; then
                sed 's/"//g' {params.output_dir}/het_fail_ind.txt | awk '{print $1, $2}' > {output.tempDir}/het_fail.txt
                plink2 --pfile {output.tempDir}/step3 --remove {output.tempDir}/het_fail.txt --make-pgen --out {output.tempDir}/step4 --threads {threads}
            else
                cp {output.tempDir}/step3.pgen {output.tempDir}/step4.pgen
                cp {output.tempDir}/step3.pvar {output.tempDir}/step4.pvar
                cp {output.tempDir}/step3.psam {output.tempDir}/step4.psam
            fi

            plink2 --pfile {output.tempDir}/step4 --make-pgen --out {output.pgen} --threads {threads}
            """

else:
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
            scripts_dir=SCRIPTS_DIR,
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
            bash {params.scripts_dir}/filterStandard.sh {output.tempDir}/pastSex {params.output_dir} {threads}
            """
