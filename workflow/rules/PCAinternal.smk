INTERNAL_PCA_METHOD = config.get("internalPCA", {}).get("method", "plink2")

rule runPcairInternalPca:
    log:
        OUT_DIR / "logs" / "runPcairInternalPca_{subset}.log",
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
    conda:
        "../../envs/ancNreport.yml"
    threads: 8
    resources:
        nodes=1,
        mem_mb=32000,
        runtime=2880,
    input:
        pgen=OUT_DIR / "{subset}" / "standardFilter.LDpruned.pgen",
        pvar=OUT_DIR / "{subset}" / "standardFilter.LDpruned.pvar",
        psam=OUT_DIR / "{subset}" / "standardFilter.LDpruned.psam",
    output:
        eigenvec=OUT_DIR / "{subset}" / "internal_pca.eigenvec",
        eigenval=OUT_DIR / "{subset}" / "internal_pca.eigenval",
        gds=temp(OUT_DIR / "{subset}" / "intermediates" / "LDpruned.gds"),
        seq_gds=temp(OUT_DIR / "{subset}" / "intermediates" / "LDpruned_seq.gds"),
        pcaobj=OUT_DIR / "{subset}" / "pcair_pcaobj.RDS",
        pcrelate=OUT_DIR / "{subset}" / "pcrelate_kinship.RDS",
        unrels=OUT_DIR / "{subset}" / "pcair_unrelated_ids.txt",
        coords=OUT_DIR / "{subset}" / "pcair_coordinates.tsv",
        grm=OUT_DIR / "{subset}" / "pcair.grm.bin",
        grmid=OUT_DIR / "{subset}" / "pcair.grm.id",
        grmN=OUT_DIR / "{subset}" / "pcair.grm.N.bin",
        plot=report(
            OUT_DIR / "{subset}" / "figures" / "pcair_pcs.svg",
            caption="../../report/pcair.rst",
            category="Internal PCA",
        ),
    params:
        out_dir=OUT_DIR / "{subset}",
        input_prefix=lambda wildcards, input: str(input.pgen)[:-4],
        color_col=config.get("internalPCA", {}).get("color_by", "None"),
        pheno_file=config.get("internalPCA", {}).get("phenotype_file", "None"),
        scripts_dir=SCRIPTS_DIR,
    shell:
        """
        if [[ "{INTERNAL_PCA_METHOD}" != "pcair" && "{INTERNAL_PCA_METHOD}" != "both" ]]; then
            echo "Skipping PC-AiR (method={INTERNAL_PCA_METHOD})"
            exit 0
        fi
        
        echo "Running PC-AiR internal PCA"
        
        mkdir -p "$(dirname {output.eigenvec})"
        
        echo "Converting pgen to bed format..."
        plink2 --pfile {params.input_prefix} --make-bed --out {params.input_prefix}
        
        Rscript {params.scripts_dir}/run_pcair_pcrelate.R \
            "{params.input_prefix}" \
            "{params.out_dir}" \
            "{params.color_col}" \
            "{params.pheno_file}" \
            "{output.plot}" \
            "{output.gds}" \
            "{output.seq_gds}"
        
        echo "Creating GRM with plink2 using unrelated samples..."
        plink2 --bfile {params.input_prefix} \
            --keep {output.unrels} \
            --make-grm-bin \
            --out {params.out_dir}/pcair_grm
        
        mv {params.out_dir}/pcair_grm.grm.bin {output.grm}
        mv {params.out_dir}/pcair_grm.grm.id {output.grmid}
        mv {params.out_dir}/pcair_grm.grm.N.bin {output.grmN}
        
        echo "PC-AiR completed"
        """


rule runPlink2ApproximatePca:
    log:
        OUT_DIR / "logs" / "runPlink2ApproximatePca_{subset}.log",
    container:
        "docker://gfanz/plink2:latest"
    conda:
        "../../envs/ancNreport.yml"
    threads: 8
    resources:
        nodes=1,
        mem_mb=32000,
        runtime=1440,
    input:
        pgen=OUT_DIR / "{subset}" / "unrelated.pgen",
        pvar=OUT_DIR / "{subset}" / "unrelated.pvar",
        psam=OUT_DIR / "{subset}" / "unrelated.psam",
    output:
        eigenvec=OUT_DIR / "{subset}" / "internal_pca_plink2.eigenvec",
        eigenval=OUT_DIR / "{subset}" / "internal_pca_plink2.eigenval",
    params:
        input_prefix=lambda wildcards, input: str(input.pgen)[:-5],
        npc=config.get("internalPCA", {}).get("npc", 20),
        tmpdir=temp(directory(OUT_DIR / "{subset}" / "intermediates" / "plink2_pca_tmp")),
    shell:
        """
        if [[ "{INTERNAL_PCA_METHOD}" != "plink2" && "{INTERNAL_PCA_METHOD}" != "both" ]]; then
            echo "Skipping PLINK2 PCA (method={INTERNAL_PCA_METHOD})"
            exit 0
        fi
        
        echo "Running PLINK2 approx PCA on unrelated samples"
        
        mkdir -p "$(dirname {output.eigenvec})"
        mkdir -p {params.tmpdir}
        
        plink2 --pfile {params.input_prefix} --make-bed --out {params.tmpdir}/unrelated
        
        plink2 --bfile {params.tmpdir}/unrelated \
            --pca approx {params.npc} \
            --out {output.eigenvec}
        
        cp {output.eigenvec}.eigenval {output.eigenval}
        
        echo "PLINK2 approx PCA completed"
        """
