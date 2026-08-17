rule crossmap:
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
    conda:
        "../../envs/ancNreport.yml"
    threads: 8
    resources:
        nodes=1,
        mem_mb=32000,
        runtime=120,
    input:
        bed=OUT_DIR / "01-Initialfilter" / "initialFilter.bed",
        bim=OUT_DIR / "01-Initialfilter" / "initialFilter.bim",
        fam=OUT_DIR / "01-Initialfilter" / "initialFilter.fam",
    output:
        # List all files that PLINK will actually create
        eigen=OUT_DIR / "04-globalAncestry" / "ref.eigenvec",
        projected=OUT_DIR / "04-globalAncestry" / "sampleRefPCscores.sscore",
        tempDir=temp(directory(OUT_DIR / "04-globalAncestry" / "intermediates")),
    params:
        method=config.get("relatedness", {}).get("method", "king"),
        grm=config.get("relatedness", {}).get("method", "king"),
        out_dir=OUT_DIR / "04-globalAncestry",
        input_prefix=OUT_DIR / "01-Initialfilter" / "initialFilter",
        input_dir=OUT_DIR / "01-Initialfilter",
        ref=REF / "1000G_highcoverage" / "1000G_highCoveragephased",
    shell:
        """

    # Since plink denote X chromosome's pseudo-autosomal region as a separate 'XY' chromosome, we want to merge to pass ontto LiftOver/CrossMap. 
    # We also reformat the numeric chromsome {1-26} to {1-22, X, Y, MT} for LiftOver/CrossMap
    plink --bfile $FILE/$NAME --merge-x no-fail --make-bed --out prep1
    plink --bfile prep1 --recode --output-chr 'MT' --out prep2
    
    rm prep.bed updated.snp updated.position updated.chr
    awk '{print $1, $4-1, $4, $2}' prep2.map > prep.bed
    
    python ${REF}/CrossMap/CrossMap.py bed ${REF}/CrossMap/GRCh37_to_GRCh38.chain.gz prep.bed study.${NAME}.lifted.bed3
    
    awk '{print $4}' study.$NAME.lifted.bed3 > updated.snp
    awk '{print $4, $3}' study.$NAME.lifted.bed3 > updated.position
    awk '{print $4, $1}' study.$NAME.lifted.bed3 > updated.chr
    plink --file prep2 --extract updated.snp --make-bed --out result1
    plink --bfile result1 --update-map updated.position --make-bed --out result2
    plink --bfile result2 --update-chr updated.chr --make-bed --out result3
    plink --bfile result3 --recode --out study.$NAME.lifted
    plink --bfile result3 --recode --make-bed --out study.$NAME.lifted
    """
