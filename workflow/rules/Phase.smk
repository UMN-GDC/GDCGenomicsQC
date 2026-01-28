
CHROMOSOMES = range(1, 22)

SLURM_LOGS = "--job-name=%x --output=logs/%x_%j.out --error=logs/%x_%j.err"
rule RFMIX:
    container: "images/my_tool.sif"  # Works for .sif, .img, or docker://
    threads: 4
    resources:
        nodes = 1,
        mem_mb = 64000,
        runtime = 1320
    input:
        bed = f"{config['OUT_DIR']}/full/unrelated.bed",
        bim = f"{config['OUT_DIR']}/full/unrelated.bim",
        fam = f"{config['OUT_DIR']}/full/unrelated.fam",
    output:
        # List all files that PLINK will actually create
        expand(f"{config['OUT_DIR']}/phased/chr{chrom}.phased.vcf.gz", chrom=CHROMOSOMES)
        expand(f"{config['OUT_DIR']}/ancestry_chr{wildcards.chrom}", chrom=CHROMOSOMES)
    params:
        method = config['relatedness']["method"],
        out_dir = f"{config['OUT_DIR']}/PCA",
        input_prefix = f"{config['OUT_DIR']}/relatedness/unrelated",
        ref= config["REF"]
    shell: """
    echo "Shapeit Phasing"

    plink --bfile $WORK/$DATATYPE/QC8 --chr $CHR --recode vcf --out ${NAME}.chr${CHR}
    bgzip -c ${NAME}.chr${CHR}.vcf > ${NAME}.chr${CHR}.vcf.gz
    bcftools index -f ${NAME}.chr${CHR}.vcf.gz
    {params.ref}/ancestry_OG/shapeit4/bin/shapeit4.2 \
        --input chr${CHR}.vcf.gz \
        --map ${REF}/ancestry_OG/chr${CHR}.b38.gmap.gz \
        --region ${CHR} \
        --output ${NAME}.chr${CHR}.phased.vcf.gz \
        --thread 8
      

    echo "RFMIX Ancestry Estimation"
    rfmix \
        -f WORK/phased/chr{wildcards.chrom}.phased.vcf.gz \
        -r {params.ref}/rfmix_ref/ALL_phase3_shapeit2_mvncall_integrated_v3plus_nounphased_rsID_genotypes_GRCh38_dbSNP.vcf.gz \
        -m {params.ref}/rfmix_ref/super_population_map_file.txt \
        -g {params.ref}/rfmix_ref/genetic_map_hg38.txt \
        -o ancestry_chr{wildcards.chrom} \
        --chromosome={wildcards.chrom}
    """
