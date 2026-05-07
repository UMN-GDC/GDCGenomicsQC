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
        pgen=OUT_DIR / "{subset}" / "initialFilter.LDpruned.pgen",
        pvar=OUT_DIR / "{subset}" / "initialFilter.LDpruned.pvar",
        psam=OUT_DIR / "{subset}" / "initialFilter.LDpruned.psam",
    output:
        pgen=OUT_DIR / "{subset}" / "unrelated.pgen",
        pvar=OUT_DIR / "{subset}" / "unrelated.pvar",
        psam=OUT_DIR / "{subset}" / "unrelated.psam",
        grm=OUT_DIR / "{subset}" / "unrelated.grm.bin",
        grmid=OUT_DIR / "{subset}" / "unrelated.grm.id",
        grmN=OUT_DIR / "{subset}" / "unrelated.grm.N.bin",
    params:
        king_cutoff=config.get("relatedness", {}).get("king_cutoff", 0.0884),
        method=config.get("relatedness", {}).get("method", "king"),
        scripts_dir=SCRIPTS_DIR,
        input_prefix=lambda wildcards, input: str(input.pgen)[:-5],
    shell:
        """

    echo "Estimating genetic relatedness"
    echo "Method: {params.method}"

    if [[ "{params.method}" == "king" || "{params.method}" == "1" ]]; then
        echo "KING ESTIMATION using PLINK2 with cutoff {params.king_cutoff}"
        plink2 --pfile {params.input_prefix} \
            --make-grm-bin \
            --king-cutoff {params.king_cutoff} \
            --king-table \
            --out {output.pgen[:-5]}_grm
        mv {output.pgen[:-5]}_grm.grm.bin {output.grm}
        mv {output.pgen[:-5]}_grm.grm.id {output.grmid}
        mv {output.pgen[:-5]}_grm.grm.N.bin {output.grmN}
        awk -v cutoff={params.king_cutoff} '$5 < cutoff {{print $1, $2}}' {output.pgen[:-5]}_grm.king.cutoff.in.id | head -n -1 > {output.pgen[:-5]}_unrel_ids.txt
        plink2 --pfile {params.input_prefix} --keep {output.pgen[:-5]}_unrel_ids.txt --make-pgen --out {output.pgen[:-5]}

    elif [[ "{params.method}" == "primus" || "{params.method}" == "2" ]]; then
        echo "PRIMUS ESTIMATION"
        mkdir -p {output.pgen[:-5]}_primus_tmp
        bash {params.scripts_dir}/run_primus.sh {params.input_prefix} {output.pgen[:-5]}_primus_tmp {config.get("REF", "/path/to/ref")}
        plink2 --bfile {output.pgen[:-5]}_primus_tmp/unrelated --make-grm-bin --out {output.pgen[:-5]}_grm
        mv {output.pgen[:-5]}_grm.grm.bin {output.grm}
        mv {output.pgen[:-5]}_grm.grm.id {output.grmid}
        mv {output.pgen[:-5]}_grm.grm.N.bin {output.grmN}
        plink2 --bfile {output.pgen[:-5]}_primus_tmp/unrelated --make-pgen --out {output.pgen[:-5]}
        rm -rf {output.pgen[:-5]}_primus_tmp

    else
        echo "ASSUMING UNRELATED: no relatedness method specified"
        cp {input.pgen} {output.pgen}
        cp {input.pvar} {output.pvar}
        cp {input.psam} {output.psam}
        plink2 --pfile {params.input_prefix} --make-grm-bin --out {output.pgen[:-5]}_grm
        mv {output.pgen[:-5]}_grm.grm.bin {output.grm}
        mv {output.pgen[:-5]}_grm.grm.id {output.grmid}
        mv {output.pgen[:-5]}_grm.grm.N.bin {output.grmN}
    fi

    """
