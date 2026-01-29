rule RFMIX:
    container: "images/my_tool.sif"  # Works for .sif, .img, or docker://
    conda: "../../envs/qcEnvironment.yml"  # Works for .sif, .img, or docker://
    threads: 8
    resources:
        nodes = 1,
        mem_mb = 64000,
        runtime = 1320
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
        ref= config["REF"]
    shell: """
    echo "Shapeit Phasing"

    plink --bfile {params.input_prefix} --chr {wildcards.CHR} --recode vcf --out {params.out_dir}/chr{wildcards.CHR}
    bgzip -c {params.out_dir}/chr{wildcards.CHR}.vcf > {params.out_dir}/chr{wildcards.CHR}.vcf.gz
    bcftools index -f {params.out_dir}/chr{wildcards.CHR}.vcf.gz
    {params.ref}/ancestry_OG/shapeit4/bin/shapeit4.2 \
        --input {params.out_dir}/chr{wildcards.CHR}.vcf.gz \
        --map {params.ref}/ancestry_OG/chr{wildcards.CHR}.b38.gmap.gz \
        --region {wildcards.CHR} \
        --output {params.out_dir}/chr{wildcards.CHR}.phased.vcf.gz \
        --thread 8
      

    echo "RFMIX Ancestry Estimation"
    rfmix \
        -f {params.out_dir}/chr{wildcards.CHR}.phased.vcf.gz \
        -r {params.ref}/rfmix_ref/ALL_phase3_shapeit2_mvncall_integrated_v3plus_nounphased_rsID_genotypes_GRCh38_dbSNP.vcf.gz \
        -m {params.ref}/rfmix_ref/super_population_map_file.txt \
        -g {params.ref}/rfmix_ref/genetic_map_hg38.txt \
        -o ancestry_chr{wildcards.CHR} \
        --chromosome={wildcards.CHR}
    """
