rule popVAE:
    conda:
        "../../envs/ancNreport.yml"
    envmodules: use("plink_module")
    threads: 8
    resources:
        nodes=1,
        mem_mb=64000,
        runtime=2880,
    input:
        pgen=OUT_DIR / "full" / "f1.b38.f2.ldpruned.pgen",
        pvar=OUT_DIR / "full" / "f1.b38.f2.ldpruned.pvar",
        psam=OUT_DIR / "full" / "f1.b38.f2.ldpruned.psam",
    output:
        f"{config.get('OUT_DIR', '/path/to/out')}/04-globalAncestry/popvae",
    params:
        input_prefix=OUT_DIR / "full" / "f1.b38.f2.ldpruned",
        npc=10,
        neighbors=50,
        ncoords=2,
        scripts_dir=SCRIPTS_DIR,
    shell:
        """
    echo "Running PopVAE"
    
    plink2 --pfile {params.input_prefix} --make-bed --out {params.input_prefix}
    plink2 --bfile {params.input_prefix} --recode vcf-iid --out {params.input_prefix}
    
    python /projects/standard/gdc/public/popvae/popvae.py \
        --infile {params.input_prefix}.vcf \
        --out {output} \
        --max_epochs 500
    """
