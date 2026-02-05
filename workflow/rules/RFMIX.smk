rule RFMIX:
    conda: "../../envs/rfmix.yml"
    threads: 8
    resources:
        nodes = 1,
        mem_mb = 64000,
        runtime = 1320,
        slurm_extra = "'--job-name=RFMIX_{wildcards.CHR}'"
    input:
        vcf = os.path.join(config['OUT_DIR'], "03-localAncestry/chr{CHR}.phased.vcf.gz"),
    output:
        # List all files that PLINK will actually create
        vcf = os.path.join(config['OUT_DIR'], "03-localAncestry/chr{CHR}_ancestry"),
    params:
        out_dir = f"{config['OUT_DIR']}/03-localAncestry",
        ref= config["REF"],
        test=config["rfmix_test"]
    shell: """
    # get common snps 
    #bcftools view -H {params.out_dir}/chr{wildcards.CHR}.phased.vcf.gz | cut -f 1,2 > {params.out_dir}/chr{wildcards.CHR}snps.txt
    # filter reference to match query
    #bcftools view -T {params.out_dir}/chr{wildcards.CHR}snps.txt ALL_phase3_shapeit2_mvncall_integrated_v3plus_nounphased_rsID_genotypes_GRCh38_dbSNP.vcf.gz -Oz -o {params.out_dir}/chr{wildcards.CHR}reference_filtered.vcf.gz
    if [ "{params.test}" = "True" ] ;  then
      echo "RFMIX Ancestry Estimation"
      rfmix \
          -f {params.out_dir}/chr{wildcards.CHR}.phased.vcf.gz \
          -r {params.ref}/rfmix_ref/ALLp3s2rsidGR38Filtered_chr{wildcards.CHR}.vcf.gz \
          -m {params.ref}/rfmix_ref/super_population_map_file.txt \
          -g {params.ref}/rfmix_ref/genetic_map_hg38.txt \
          -o chr{wildcards.CHR}_ancestry \
          -e 1 \
          -t 2 \
          --n-threads {threads} \
          --chromosome={wildcards.CHR}
    else
      rfmix \
          -f {params.out_dir}/chr{wildcards.CHR}.phased.vcf.gz \
          -r {params.ref}/rfmix_ref/ALLp3s2rsidGR38Filtered_chr{wildcards.CHR}.vcf.gz \
          -m {params.ref}/rfmix_ref/super_population_map_file.txt \
          -g {params.ref}/rfmix_ref/genetic_map_hg38.txt \
          -o chr{wildcards.CHR}_ancestry \
          --n-threads {threads} \
          --chromosome={wildcards.CHR}
    fi
    """
