rule PCAinternal:
    container: "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
    conda: "../../envs/ancNreport.yml"
    threads: 8
    resources:
        nodes = 1,
        mem_mb = 32000,
        runtime = 2880,
    input:
        pgen = OUT_DIR / "{subset}" / "standardFilter.LDpruned.pgen",
        pvar = OUT_DIR / "{subset}" / "standardFilter.LDpruned.pvar",
        psam = OUT_DIR / "{subset}" / "standardFilter.LDpruned.psam",
    output:
        gds = temp(OUT_DIR / "{subset}" / "intermediates" / "LDpruned.gds"),
        seq_gds = temp(OUT_DIR / "{subset}" / "intermediates" / "LDpruned_seq.gds"),
        pcaobj = OUT_DIR / "{subset}" / "pcair_pcaobj.RDS",
        unrels = OUT_DIR / "{subset}" / "pcair_unrelated_ids.txt",
        pcrelate = OUT_DIR / "{subset}" / "pcrelate_kinship.RDS",
        coords = OUT_DIR / "{subset}" / "pcair_coordinates.tsv",
        plot = report(OUT_DIR / "{subset}" / "figures" / "pcair_pcs.svg", 
                      caption = "../../report/pcair.rst", 
                      category = "Internal PCA"),
    params:
        out_dir = OUT_DIR / "{subset}",
        input_prefix = lambda wildcards, input: str(input.pgen)[:-4],
        color_col = config.get('internalPCA', {}).get('color_by', 'None'),
        pheno_file = config.get('internalPCA', {}).get('phenotype_file', 'None'),
    shell: 
        """
        echo "Running PC-AiR and PC-Relate for internal samples"
        
        mkdir -p "$(dirname {output.plot})"
        
        echo "Converting pgen to bed format..."
        plink2 --pfile {params.input_prefix} --export bed --out {params.input_prefix}
        
        Rscript scripts/run_pcair_pcrelate.R \
            "{params.input_prefix}" \
            "{params.out_dir}" \
            "{params.color_col}" \
            "{params.pheno_file}" \
            "{output.plot}" \
            "{output.gds}" \
            "{output.seq_gds}"
        
        echo "PC-AiR and PC-Relate completed successfully"
        """
