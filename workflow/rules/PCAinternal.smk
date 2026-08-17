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
        unrels=OUT_DIR / "{subset}" / "pcair_unrelated_ids.txt",
        coords=OUT_DIR / "{subset}" / "pcair_coordinates.tsv",
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
        echo "Running PC-AiR internal PCA"
        
        mkdir -p "$(dirname {output.eigenvec})"
        
        echo "Converting pgen to bed format..."
        plink2 --pfile {params.input_prefix} --export bed --out {params.input_prefix}
        
        Rscript {params.scripts_dir}/run_pcair_pcrelate.R \
            "{params.input_prefix}" \
            "{params.out_dir}" \
            "{params.color_col}" \
            "{params.pheno_file}" \
            "{output.plot}" \
            "{output.gds}" \
            "{output.seq_gds}"
        
        echo "PC-AiR completed"
        """
rule runPlink2ApproximatePca:
    log:
        OUT_DIR / "logs" / "runPlink2ApproximatePca_{subset}.log",
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
    conda:
        "../../envs/ancNreport.yml"
    threads: 8
    resources:
        nodes=1,
        mem_mb=32000,
        runtime=1440,
    input:
        pgen=OUT_DIR / "{subset}" / "standardFilter.LDpruned.pgen",
        pvar=OUT_DIR / "{subset}" / "standardFilter.LDpruned.pvar",
        psam=OUT_DIR / "{subset}" / "standardFilter.LDpruned.psam",
        keep=OUT_DIR / "{subset}" / "unrelated.keep",
    output:
        eigenvec=OUT_DIR / "{subset}" / "internal_pca_plink2.eigenvec",
        eigenval=OUT_DIR / "{subset}" / "internal_pca_plink2.eigenval",
    params:
        input_prefix=lambda wc, input: str(input.pgen)[:-5],
        output_prefix=lambda wc, output: str(output.eigenvec).replace(
            ".eigenvec", ""
        ),
        npc=config.get("internalPCA", {}).get("npc", 10),
    shell:
        """
        set -euo pipefail

        plink2 --pfile {params.input_prefix} \
          --keep {input.keep} \
          --pca approx {params.npc} \
          --out {params.output_prefix} \
          --threads {threads}
        """
