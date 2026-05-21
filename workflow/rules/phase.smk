def get_chrom(wildcards):
    return wildcards.CHR.replace('.phased', '')


def get_input_pgen(wildcards):
    if "{CHR}" in config.get("INPUT", ""):
        return OUT_DIR / "full" / f"f1.f2_{wildcards.CHR}.pgen"
    else:
        return OUT_DIR / "full" / "f1.b38.f2.pgen"


def get_input_pvar(wildcards):
    if "{CHR}" in config.get("INPUT", ""):
        return OUT_DIR / "full" / f"f1.f2_{wildcards.CHR}.pvar"
    else:
        return OUT_DIR / "full" / "f1.b38.f2.pvar"


def get_input_psam(wildcards):
    if "{CHR}" in config.get("INPUT", ""):
        return OUT_DIR / "full" / f"f1.f2_{wildcards.CHR}.psam"
    else:
        return OUT_DIR / "full" / "f1.b38.f2.psam"


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
            OUT_DIR / "full" / f"f1.f2_{wildcards.CHR}.pgen"
            if "{CHR}" in config.get("INPUT", "")
            else OUT_DIR / "full" / "f1.b38.f2.pgen"
        ),
        pgen=get_input_pgen,
        pvar=get_input_pvar,
        psam=get_input_psam,
    output:
        vcf=OUT_DIR / "02-localAncestry" / "chr{CHR}.vcf.gz",
        csi=OUT_DIR / "02-localAncestry" / "chr{CHR}.vcf.gz.csi",
    params:
        out_dir=OUT_DIR / "02-localAncestry",
        input_prefix=lambda wildcards, input: input.pgen[:-5],
        chrom=get_chrom,
    shell:
        """
        plink2 --pfile {params.input_prefix} --chr {params.chrom} --allow-extra-chr --recode vcf bgz --out {params.out_dir}/chr{wildcards.CHR}
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
          plink2 --vcf {input.vcf} --bp-space 100000 --thin-indiv {params.thin} --recode vcf bgz --out {params.out_dir}/chr{wildcards.CHR}.thinned
          mv {params.out_dir}/chr{wildcards.CHR}.thinned.vcf.gz {params.out_dir}/chr{wildcards.CHR}.vcf.gz
          bcftools index -f {params.out_dir}/chr{wildcards.CHR}.vcf.gz
          echo "Running shapeit4 in test mode"
          shapeit4 \
              --input {params.out_dir}/chr{wildcards.CHR}.vcf.gz \
              --map {input.gmap} \
              --reference {input.ref} \
              --region {params.chrom} \
              --log {params.out_dir}/chr{wildcards.CHR}.phased.log \
              --thread {threads} \
              --mcmc-iterations 1b,1p,1m \
              --output {params.out_dir}/chr{wildcards.CHR}.phased.vcf
        else
          shapeit4 \
              --input {input.vcf} \
              --map {input.gmap} \
              --reference {input.ref} \
              --region {params.chrom} \
              --log {params.out_dir}/chr{wildcards.CHR}.phased.log \
              --thread {threads} \
              --output {params.out_dir}/chr{wildcards.CHR}.phased.vcf
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
