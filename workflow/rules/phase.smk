def get_chrom(wildcards):
    return wildcards.CHR.replace('.phased', '')


def get_input_pgen(wildcards):
    if "{CHR}" in config.get("INPUT", ""):
        return OUT_DIR / "full" / f"initialFilter_{get_chrom(wildcards)}.pgen"
    else:
        return OUT_DIR / "full" / "initialFilter.pgen"


def get_input_pvar(wildcards):
    if "{CHR}" in config.get("INPUT", ""):
        return OUT_DIR / "full" / f"initialFilter_{get_chrom(wildcards)}.pvar"
    else:
        return OUT_DIR / "full" / "initialFilter.pvar"


def get_input_psam(wildcards):
    if "{CHR}" in config.get("INPUT", ""):
        return OUT_DIR / "full" / f"initialFilter_{get_chrom(wildcards)}.psam"
    else:
        return OUT_DIR / "full" / "initialFilter.psam"


rule convertPgenToVcf:
    log:
        OUT_DIR / "logs" / "Convert_{CHR}.log",
    container: "oras://ghcr.io/coffm049/gdcgenomicsqc/rfmix:latest"
    conda: "../../envs/rfmix.yml"
    threads: 8
    resources:
        nodes=1,
        mem_mb=16000,
        runtime=120,
    input:
        std=OUT_DIR / "full" / "standardFilter.pgen",
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
        plink2 --pfile {params.input_prefix} --chr {params.chrom} --recode vcf bgz --out {params.out_dir}/chr{wildcards.CHR}
        bcftools index -f {params.out_dir}/chr{wildcards.CHR}.vcf.gz
        """


rule phaseWithShapeit:
    log:
        OUT_DIR / "logs" / "Phase_{CHR}.log",
    container: "oras://ghcr.io/coffm049/gdcgenomicsqc/rfmix:latest"
    conda: "../../envs/rfmix.yml"
    threads: 8
    resources:
        nodes=1,
        mem_mb=64000,
        runtime=1320,
    input:
        vcf=OUT_DIR / "02-localAncestry" / "chr{CHR}.vcf.gz",
        gmap=ancient(REF / "gmaps" / "chr{CHR}.b38.gmap.gz"),
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
              --region {params.chrom} \
              --log {params.out_dir}/chr{wildcards.CHR}.phased.log \
              --thread {threads} \
              --mcmc-iterations 1b,1p,1m \
              --output {params.out_dir}/chr{wildcards.CHR}.phased.vcf
        else
          shapeit4 \
              --input {input.vcf} \
              --map {input.gmap} \
              --region {params.chrom} \
              --log {params.out_dir}/chr{wildcards.CHR}.phased.log \
              --thread {threads} \
              --output {params.out_dir}/chr{wildcards.CHR}.phased.vcf
        fi
        """


rule compressAndIndexVcf:
    log:
        OUT_DIR / "logs" / "Compress_{CHR}.log",
    container: "oras://ghcr.io/coffm049/gdcgenomicsqc/rfmix:latest"
    conda: "../../envs/rfmix.yml"
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
