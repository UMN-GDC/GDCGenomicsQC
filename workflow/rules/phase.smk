rule Phase:
    container: "oras://ghcr.io/coffm049/gdcgenomicsqc/rfmix:latest"
    conda: "../../envs/rfmix.yml"
    threads: 8
    resources:
        nodes = 1,
        mem_mb = 64000,
        runtime = 1320,
    input:
        pgen =   OUT_DIR / "full" / "initialFilter.pgen",
        pvar =   OUT_DIR / "full" / "initialFilter.pvar",
        psam =   OUT_DIR / "full" / "initialFilter.psam",
    output:
        # List all files that PLINK will actually create
        vcf = OUT_DIR / "02-localAncestry" / "chr{CHR}.phased.vcf.gz",
    params:
        out_dir = OUT_DIR / "02-localAncestry",
        input_prefix = lambda wildcards, input: input.pgen[:-5],
        ref= config["REF"],
        test = config["localAncestry"]["test"],
        thin = config["localAncestry"]["thin_subjects"],
    shell: """
    echo "Shapeit Phasing"

    plink2 --pfile {params.input_prefix} --chr {wildcards.CHR} --recode vcf bgz --out {params.out_dir}/chr{wildcards.CHR}
    bcftools index -f {params.out_dir}/chr{wildcards.CHR}.vcf.gz
    
    if [ "{params.test}" = "True" ] ;  then
      plink2 --vcf {params.out_dir}/chr{wildcards.CHR}.vcf.gz --bp-space 100000 --thin-indiv {params.thin} --recode vcf bgz --out {params.out_dir}/chr{wildcards.CHR}.thinned
      mv {params.out_dir}/chr{wildcards.CHR}.thinned.vcf.gz {params.out_dir}/chr{wildcards.CHR}.vcf.gz
      bcftools index -f {params.out_dir}/chr{wildcards.CHR}.vcf.gz
      echo "{params.out_dir}/chr{wildcards.CHR}.phased.vcf.gz"
      shapeit4 \
          --input {params.out_dir}/chr{wildcards.CHR}.vcf.gz \
          --map {params.ref}/ancestry_OG/chr{wildcards.CHR}.b38.gmap.gz \
          --region {wildcards.CHR} \
          --log {params.out_dir}/chr{wildcards.CHR}.phased.log \
          --thread {threads} \
          --mcmc-iterations 1b,1p,1m \
          --output {params.out_dir}/chr{wildcards.CHR}.phased.vcf
      bgzip -c {params.out_dir}/chr{wildcards.CHR}.phased.vcf > {params.out_dir}/chr{wildcards.CHR}.phased.vcf.gz
      bcftools index -f {params.out_dir}/chr{wildcards.CHR}.phased.vcf.gz
      rm {params.out_dir}/chr{wildcards.CHR}.phased.vcf
      
    else
      shapeit4 \
          --input {params.out_dir}/chr{wildcards.CHR}.vcf.gz \
          --map {params.ref}/ancestry_OG/chr{wildcards.CHR}.b38.gmap.gz \
          --region {wildcards.CHR} \
          --log {params.out_dir}/chr{wildcards.CHR}.phased.log \
          --thread {threads} \
          --output {params.out_dir}/chr{wildcards.CHR}.phased.vcf
      bgzip -c {params.out_dir}/chr{wildcards.CHR}.phased.vcf > {params.out_dir}/chr{wildcards.CHR}.phased.vcf.gz
      bcftools index -f {params.out_dir}/chr{wildcards.CHR}.phased.vcf.gz
      rm {params.out_dir}/chr{wildcards.CHR}.phased.vcf
    fi
      

    """
