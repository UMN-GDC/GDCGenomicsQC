rule checkRelatednessExtractUnrelated:
    log:
        OUT_DIR / "logs" / "checkRelatednessExtractUnrelated_{subset}.log",
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
    conda:
        "../../envs/ancNreport.yml"
    resources:
        nodes=1,
        mem_mb=128000,
        runtime=720,
    input:
        bed=OUT_DIR / "{subset}" / "initialFilter.LDpruned.bed",
        bim=OUT_DIR / "{subset}" / "initialFilter.LDpruned.bim",
        fam=OUT_DIR / "{subset}" / "initialFilter.LDpruned.fam",
    output:
        bed=OUT_DIR / "{subset}" / "unrelated.bed",
        bim=OUT_DIR / "{subset}" / "unrelated.bim",
        fam=OUT_DIR / "{subset}" / "unrelated.fam",
        grm=OUT_DIR / "{subset}" / "unrelated.grm.bin",
        grmid=OUT_DIR / "{subset}" / "unrelated.grm.id",
        grmN=OUT_DIR / "{subset}" / "unrelated.grm.N.bin",
    params:
        king_cutoff=config.get("relatedness", {}).get("king_cutoff", 0.0884),
        method=config.get("relatedness", {}).get("method", "king"),
    shell:
        """
    
    echo "Estimating genetic relatedness"
    echo "Method: {params.method}"

    if [[ "{params.method}" == "king" || "{params.method}" == "1" ]]; then
        echo "KING ESTIMATION using PLINK2 with cutoff {params.king_cutoff}"
        plink2 --bfile {input.bed[:-4]} \
            --make-grm-bin \
            --king-cutoff {params.king_cutoff} \
            --make-bed \
            --out {output.bed[:-4]}
        mv {output.bed[:-4]}.grm.bin {output.grm}
        mv {output.bed[:-4]}.grm.id {output.grmid}
        mv {output.bed[:-4]}.grm.N.bin {output.grmN}
    elif [[ "{params.method}" == "primus" || "{params.method}" == "2" ]]; then
        echo "PRIMUS ESTIMATION"
        bash scripts/run_primus.sh {input.bed[:-4]} {output.bed[:-4]}
        plink2 --bfile {output.bed[:-4]} --make-grm-bin --out {output.bed[:-4]}
        mv {output.bed[:-4]}.grm.bin {output.grm}
        mv {output.bed[:-4]}.grm.id {output.grmid}
        mv {output.bed[:-4]}.grm.N.bin {output.grmN}
    else
        echo "ASSUMING UNRELATED: no relatedness method specified"
        cp {input.bed} {output.bed}
        cp {input.bim} {output.bim}
        cp {input.fam} {output.fam}
        plink2 --bfile {input.bed[:-4]} --make-grm-bin --out {output.bed[:-4]}
        mv {output.bed[:-4]}.grm.bin {output.grm}
        mv {output.bed[:-4]}.grm.id {output.grmid}
        mv {output.bed[:-4]}.grm.N.bin {output.grmN}
    fi

    """
