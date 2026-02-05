rule Phase:
    conda: "../../envs/rfmix.yml"
    threads: 8
    resources:
        nodes = 1,
        mem_mb = 64000,
        runtime = 1320,
    input:
        bed =   os.path.join(config['OUT_DIR'], "02-relatedness/standardFiltered.bed"),
        bim =   os.path.join(config['OUT_DIR'], "02-relatedness/standardFiltered.bim"),
        fam =   os.path.join(config['OUT_DIR'], "02-relatedness/standardFiltered.fam"),
    output:
        # List all files that PLINK will actually create
        vcf = os.path.join(config['OUT_DIR'], "03-localAncestry/chr{CHR}.phased.vcf.gz"),
    params:
        out_dir = f"{config['OUT_DIR']}/03-localAncestry",
        input_prefix = lambda wildcards, input: input.bed[:-4],
        ref= config["REF"],
        test=config["rfmix_test"]
    shell: """
    echo "Shapeit Phasing"

    plink2 --bfile {params.input_prefix} --chr {wildcards.CHR} --recode vcf bgz --out {params.out_dir}/chr{wildcards.CHR}
    bcftools index -f {params.out_dir}/chr{wildcards.CHR}.vcf.gz
    
    module load shapeit/4.2.2
    if [ "{params.test}" = "True" ] ;  then
      echo "{params.out_dir}/chr{wildcards.CHR}.phased.vcf.gz"
      shapeit \
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
      shapeit \
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
