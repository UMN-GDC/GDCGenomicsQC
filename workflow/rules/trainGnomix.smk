rule train_gnomix:
    conda: "../../envs/gnomix.yml"
    threads: 4
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
        phase=True
    shell: """

    python3 {params.ref}/gnomix/gnomix.py \
      {params.ref}/1000G_GRCh38/ALL.chr{wildcards.CHR}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz \
      {params.out_dir} \
      {wildcards.CHR} \
      {params.phase} \
      {params.ref}/rfmix_ref/genetic_map_hg38.txt \
      {params.ref}/1000G_GRCh38/ALL.chr{wildcards.CHR}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz \
      {params.ref}/1000G_GRCh38/1000G.popu \
      {params.ref}/gnomix/config.yaml

    """

