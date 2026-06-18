def get_chrom(wildcards):
    return wildcards.CHR.replace('.phased', '')


def get_input_pgen(wildcards):
    if "{CHR}" in config.get("INPUT", ""):
        return OUT_DIR / "full" / f"f1.f2_{wildcards.CHR}.unrel.pgen"
    else:
        return OUT_DIR / "full" / "f1.b38.f2.unrel.pgen"


def get_input_pvar(wildcards):
    if "{CHR}" in config.get("INPUT", ""):
        return OUT_DIR / "full" / f"f1.f2_{wildcards.CHR}.unrel.pvar"
    else:
        return OUT_DIR / "full" / "f1.b38.f2.unrel.pvar"


def get_input_psam(wildcards):
    if "{CHR}" in config.get("INPUT", ""):
        return OUT_DIR / "full" / f"f1.f2_{wildcards.CHR}.unrel.psam"
    else:
        return OUT_DIR / "full" / "f1.b38.f2.unrel.psam"


if INPUT_IS_PER_CHROMOSOME:
    rule extractUnrelatedSubjects:
        log:
            OUT_DIR / "logs" / "extractUnrelated_{CHR}.log",
        container:
            "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
        conda:
            "../../envs/ancNreport.yml"
        envmodules: *([config.get("plink_module")] if config.get("plink_module") else [])
        threads: 4
        resources:
            nodes=1,
            mem_mb=16000,
            runtime=60,
        input:
            pgen=OUT_DIR / "full" / "f1.f2_{CHR}.pgen",
            pvar=OUT_DIR / "full" / "f1.f2_{CHR}.pvar",
            psam=OUT_DIR / "full" / "f1.f2_{CHR}.psam",
            keep=OUT_DIR / "full" / "f1.b38.ldpruned.unrelated_ids.txt",
        output:
            pgen=OUT_DIR / "full" / "f1.f2_{CHR}.unrel.pgen",
            pvar=OUT_DIR / "full" / "f1.f2_{CHR}.unrel.pvar",
            psam=OUT_DIR / "full" / "f1.f2_{CHR}.unrel.psam",
        params:
            input_prefix=lambda wildcards, input: input.pgen[:-5],
            output_prefix=lambda wildcards, output: output.pgen[:-5],
        shell:
            """
            plink2 --pfile {params.input_prefix} --keep {input.keep} --make-pgen --out {params.output_prefix} --threads {threads}
            """

else:
    rule extractUnrelatedSubjects:
        log:
            OUT_DIR / "logs" / "extractUnrelated.log",
        container:
            "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
        conda:
            "../../envs/ancNreport.yml"
        envmodules: *([config.get("plink_module")] if config.get("plink_module") else [])
        threads: 4
        resources:
            nodes=1,
            mem_mb=16000,
            runtime=60,
        input:
            pgen=OUT_DIR / "full" / "f1.b38.f2.pgen",
            pvar=OUT_DIR / "full" / "f1.b38.f2.pvar",
            psam=OUT_DIR / "full" / "f1.b38.f2.psam",
            keep=OUT_DIR / "full" / "f1.b38.ldpruned.unrelated_ids.txt",
        output:
            pgen=OUT_DIR / "full" / "f1.b38.f2.unrel.pgen",
            pvar=OUT_DIR / "full" / "f1.b38.f2.unrel.pvar",
            psam=OUT_DIR / "full" / "f1.b38.f2.unrel.psam",
        params:
            input_prefix=lambda wildcards, input: input.pgen[:-5],
            output_prefix=lambda wildcards, output: output.pgen[:-5],
        shell:
            """
            plink2 --pfile {params.input_prefix} --keep {input.keep} --make-pgen --out {params.output_prefix} --threads {threads}
            """


rule convertPgenToVcf:
    log:
        OUT_DIR / "logs" / "Convert_{CHR}.log",
    container: "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
    conda: "../../envs/rfmix.yml"
    envmodules: *[m for m in (config.get("plink_module"), config.get("bcftools_module")) if m]
    threads: 8
    resources:
        nodes=1,
        mem_mb=16000,
        runtime=120,
    input:
        std=lambda wildcards: (
            OUT_DIR / "full" / f"f1.f2_{wildcards.CHR}.unrel.pgen"
            if "{CHR}" in config.get("INPUT", "")
            else OUT_DIR / "full" / "f1.b38.f2.unrel.pgen"
        ),
        pgen=get_input_pgen,
        pvar=get_input_pvar,
        psam=get_input_psam,
        ref=ancient(REF / "1000G_highcoverage" / "1kGP_high_coverage_Illumina.chr{CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"),
    output:
        vcf=OUT_DIR / "02-localAncestry" / "chr{CHR}.vcf.gz",
        csi=OUT_DIR / "02-localAncestry" / "chr{CHR}.vcf.gz.csi",
    params:
        out_dir=OUT_DIR / "02-localAncestry",
        input_prefix=lambda wildcards, input: input.pgen[:-5],
        chrom=get_chrom,
    shell:
        """
        plink2 --pfile {params.input_prefix} --chr {params.chrom} --allow-extra-chr --make-pgen --out {params.out_dir}/chr{wildcards.CHR}.temp --set-all-var-ids @:#:\\$r:\\$a --snps-only just-acgt
        awk '!/^#/ && (($4=="A" && $5=="T") || ($4=="T" && $5=="A") || ($4=="C" && $5=="G") || ($4=="G" && $5=="C")) {print $3}' {params.out_dir}/chr{wildcards.CHR}.temp.pvar > {params.out_dir}/chr{wildcards.CHR}.palindromic_snps.txt
        plink2 --pfile {params.out_dir}/chr{wildcards.CHR}.temp \
                       --exclude {params.out_dir}/chr{wildcards.CHR}.palindromic_snps.txt \
                       --output-chr chrM \
                       --export vcf bgz \
                       --out {params.out_dir}/chr{wildcards.CHR}
        rm {params.out_dir}/chr{wildcards.CHR}.temp.* {params.out_dir}/chr{wildcards.CHR}.palindromic_snps.txt
        bcftools index -f {params.out_dir}/chr{wildcards.CHR}.vcf.gz
        bcftools isec -p {params.out_dir}/chr{wildcards.CHR}.strict -n =2 {params.out_dir}/chr{wildcards.CHR}.vcf.gz {input.ref}
        mv {params.out_dir}/chr{wildcards.CHR}.strict/0002.vcf.gz {params.out_dir}/chr{wildcards.CHR}.vcf.gz
        rm -rf {params.out_dir}/chr{wildcards.CHR}.strict
        bcftools index -f {params.out_dir}/chr{wildcards.CHR}.vcf.gz
        """


rule phaseWithShapeit:
    log:
        OUT_DIR / "logs" / "Phase_{CHR}.log",
    container: "oras://ghcr.io/coffm049/gdcgenomicsqc/rfmix:v1"
    conda: "../../envs/rfmix.yml"
    envmodules: *[m for m in (config.get("plink_module"), config.get("bcftools_module"), config.get("shapeit_module")) if m]
    threads: 8
    resources:
        nodes=1,
        mem_mb=64000,
        runtime=1320,
    input:
        vcf=OUT_DIR / "02-localAncestry" / "chr{CHR}.vcf.gz",
        ref=ancient(REF / "1000G_highcoverage" / "1kGP_high_coverage_Illumina.chr{CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"),
        gmap=ancient(REF / "gmaps" / "hg38map.chr{CHR}.txt"),
    output:
        vcf=temp(OUT_DIR / "02-localAncestry" / "chr{CHR}.phased.vcf"),
    params:
        out_dir=OUT_DIR / "02-localAncestry",
        test=config.get("localAncestry", {}).get("test", False),
        thin=config.get("localAncestry", {}).get("thin_subjects", 0.1),
        chrom=get_chrom,
    shell:
        """
        echo "Shapeit Phasing"

        if [ "{params.test}" = "True" ] ; then
          plink2 --vcf {input.vcf} --bp-space 100000 --thin-indiv {params.thin} --export vcf bgz --out {params.out_dir}/chr{wildcards.CHR}.thinned
          mv {params.out_dir}/chr{wildcards.CHR}.thinned.vcf.gz {params.out_dir}/chr{wildcards.CHR}.vcf.gz
          bcftools index -f {params.out_dir}/chr{wildcards.CHR}.vcf.gz
          echo "Running shapeit4 in test mode"
          awk '{{print "chr" $0}}' {input.gmap} > {params.out_dir}/chr{wildcards.CHR}.fixed_map.txt
          shapeit4 \
              --input {params.out_dir}/chr{wildcards.CHR}.vcf.gz \
              --map {params.out_dir}/chr{wildcards.CHR}.fixed_map.txt \
              --region chr{params.chrom} \
              --log {params.out_dir}/chr{wildcards.CHR}.phased.log \
              --thread {threads} \
              --mcmc-iterations 1b,1p,1m \
              --output {output.vcf} \
              --reference {input.ref}
          rm -f {params.out_dir}/chr{wildcards.CHR}.fixed_map.txt
        else
          awk '{{print "chr" $0}}' {input.gmap} > {params.out_dir}/chr{wildcards.CHR}.fixed_map.txt
          shapeit4 \
              --input {input.vcf} \
              --map {params.out_dir}/chr{wildcards.CHR}.fixed_map.txt \
              --region chr{params.chrom} \
              --log {params.out_dir}/chr{wildcards.CHR}.phased.log \
              --thread {threads} \
              --output {output.vcf} \
              --reference {input.ref}
          rm -f {params.out_dir}/chr{wildcards.CHR}.fixed_map.txt
        fi
        """


rule compressAndIndexVcf:
    log:
        OUT_DIR / "logs" / "Compress_{CHR}.log",
    container: "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
    conda: "../../envs/rfmix.yml"
    envmodules: *([config.get("bcftools_module")] if config.get("bcftools_module") else [])
    threads: 4
    resources:
        nodes=1,
        mem_mb=16000,
        runtime=60,
    input:
        vcf=OUT_DIR / "02-localAncestry" / "chr{CHR}.phased.vcf",
    output:
        vcf=OUT_DIR / "02-localAncestry" / "chr{CHR}.phased.vcf.gz",
        csi=OUT_DIR / "02-localAncestry" / "chr{CHR}.phased.vcf.gz.csi",
    params:
        out_dir=OUT_DIR / "02-localAncestry",
        chrom=get_chrom,
    shell:
        """
        bgzip -c {input.vcf} > {params.out_dir}/chr{wildcards.CHR}.phased.vcf.gz
        bcftools index -f {params.out_dir}/chr{wildcards.CHR}.phased.vcf.gz
        rm {input.vcf}
        """
