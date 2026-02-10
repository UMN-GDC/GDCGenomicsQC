rule train_gnomix:
    conda: "../../envs/gnomix.yml"
    threads: 2
    resources:
        nodes = 1,
        mem_mb = 64000,
        runtime = 1320,
    input:
        vcf = os.path.join(config['OUT_DIR'], "03-localAncestry/chr{CHR}.phased.vcf.gz"),
    output:
        # List all files that PLINK will actually create
        vcf = os.path.join(config['OUT_DIR'], "03-localAncestry/chr{CHR}_ancestry_gnomix"),
    params:
        out_dir = f"{config['OUT_DIR']}/03-localAncestry",
        ref= config["REF"],
        test=config["rfmix_test"],
        PHASE="False"
    shell: """

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
      {params.ref}/1000G_GRCh38/ALL.chr{wildcards.CHR}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz \
      {params.out_dir} \
      {wildcards.CHR} \
      {params.PHASE} \
      {params.ref}/rfmix_ref/genetic_map_hg38.txt \
      {params.ref}/1000G_GRCh38/ALL.chr{wildcards.CHR}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz \
      {params.ref}/1000G_GRCh38/1000G.popu \
      {params.out_dir}/chr{wildcards.CHR}config.yaml
    """

