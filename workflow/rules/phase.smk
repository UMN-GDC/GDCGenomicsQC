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
        pgen=OUT_DIR / "full" / "initialFilter.pgen",
        pvar=OUT_DIR / "full" / "initialFilter.pvar",
        psam=OUT_DIR / "full" / "initialFilter.psam",
    output:
        vcf=OUT_DIR / "02-localAncestry" / "chr{CHR}.vcf.gz",
        csi=OUT_DIR / "02-localAncestry" / "chr{CHR}.vcf.gz.csi",
    params:
        out_dir=OUT_DIR / "02-localAncestry",
        input_prefix=lambda wildcards, input: input.pgen[:-5],
    shell:
        """
        plink2 --pfile {params.input_prefix} --chr {wildcards.CHR} --recode vcf bgz --out {params.out_dir}/chr{wildcards.CHR}
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
        gmap=REF / "gmaps" / "chr{CHR}.b38.gmap.gz",
    output:
        vcf=temp(OUT_DIR / "02-localAncestry" / "chr{CHR}.phased.vcf"),
    params:
        out_dir=OUT_DIR / "02-localAncestry",
        test=config.get("localAncestry", {}).get("test", False),
        thin=config.get("localAncestry", {}).get("thin_subjects", 0.1),
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
              --region {wildcards.CHR} \
              --log {params.out_dir}/chr{wildcards.CHR}.phased.log \
              --thread {threads} \
              --mcmc-iterations 1b,1p,1m \
              --output {params.out_dir}/chr{wildcards.CHR}.phased.vcf
        else
          shapeit4 \
              --input {input.vcf} \
              --map {input.gmap} \
              --region {wildcards.CHR} \
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
    shell:
        """
        bgzip -c {input.vcf} > {params.out_dir}/chr{wildcards.CHR}.phased.vcf.gz
        bcftools index -f {params.out_dir}/chr{wildcards.CHR}.phased.vcf.gz
        rm {input.vcf}
        """