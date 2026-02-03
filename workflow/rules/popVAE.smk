rule popVAE:
    conda: "../../envs/ancNreport.yml"
    threads: 8
    resources:
        nodes = 1,
        mem_mb = 64000,
        runtime = 2880,
        slurm_extra = "'--job-name=UMAP'"
    input:
        bed = os.path.join(config['OUT_DIR'], "02-relatedness/standardFiltered.LDpruned.bed"),
        bim = os.path.join(config['OUT_DIR'], "02-relatedness/standardFiltered.LDpruned.bim"),
        fam = os.path.join(config['OUT_DIR'], "02-relatedness/standardFiltered.LDpruned.fam"),
    output:
        f"{config['OUT_DIR']}/04-globalAncestry/popvae",
    params:
        input_prefix = lambda wildcards, input: input.bed[:-4],
        npc = 10,
        neighbors = 50,
        ncoords = 2,
    shell: """
    echo "Running PopVAE"
    
    plink2 --bfile {params.input_prefix}  --recode vcf-iid --out {params.input_prefix}
    
    # Works Run PopVAE using the specific Python from the popvae environment
    python /projects/standard/gdc/public/popvae/popvae.py \
        --infile {params.input_prefix}.vcf \
        --out {output} \
        --max_epochs 500
    """
