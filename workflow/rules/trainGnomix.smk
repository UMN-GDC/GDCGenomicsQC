rule train_gnomix:
    conda:
        "../../envs/gnomix.yml"
    threads: 2
    resources:
        nodes=1,
        mem_mb=128000,
        runtime=1320,
    input:
        vcf=os.path.join(config.get("OUT_DIR", "/path/to/out"), "03-localAncestry/chr{CHR}.phased.vcf.gz"),
    output:
        vcf=os.path.join(config.get("OUT_DIR", "/path/to/out"), "03-localAncestry/chr{CHR}_ancestry_gnomix"),
    params:
        out_dir=f"{config.get('OUT_DIR', '/path/to/out')}/03-localAncestry",
        ref=config.get("REF", "/path/to/ref"),
        test=config.get("localAncestry", {}).get("test", False),
        PHASE="False",
    shell:
        """

    # make a quick runs possible
    if [ "{params.test}" = "True" ] ;  then
      sed -e 's/name: model/name: chr{wildcards.CHR}_model/;' \
          -e 's/inference: /inference: fast/;' \
          -e 's/gens: [0, 2, 4, 6, 8, 12, 16, 24]/gens: [0]/;' \
          -e 's/window_size_cM: 0.2/window_size_cM: 0.1/;' \
          -e 's/smooth_size: 75/smooth_size: 5/;' \
          -e 's/context_ratio: 0.5/context_ratio: 0.1/;' \
          -e 's/retrain_base: True/retrain_base: False/;' \
          {params.ref}/gnomix/config.yaml > {params.out_dir}/chr{wildcards.CHR}config.yaml
    else 
      sed -e 's/name: model/name: chr{wildcards.CHR}_model/' {params.ref}/gnomix/config.yaml > {params.out_dir}/chr{wildcards.CHR}config.yaml
    fi
    
    python3 {params.ref}/gnomix/gnomix.py \
      {params.ref}/1000G_highcoverage/1kGP_high_coverage_Illumina.chr{wildcards.CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
      {params.out_dir} \
      {wildcards.CHR} \
      {params.PHASE} \
      {params.ref}/rfmix_ref/genetic_map_hg38.txt \
      {params.ref}/1000G_highcoverage/1kGP_high_coverage_Illumina.chr{wildcards.CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
      {params.ref}/1000G_highcoverage/population.txt \
      {params.out_dir}/chr{wildcards.CHR}config.yaml
    """
