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
        output_prefix=lambda wildcards, output: str(output.pgen)[:-5],
        ref_path=config.get("REF", "/path/to/ref"),
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
            --out {params.output_prefix}_grm
        mv {params.output_prefix}_grm.grm.bin {output.grm}
        mv {params.output_prefix}_grm.grm.id {output.grmid}
        mv {params.output_prefix}_grm.grm.N.bin {output.grmN}
        awk -v cutoff={params.king_cutoff} '$5 < cutoff {{print $1, $2}}' {params.output_prefix}_grm.king.cutoff.in.id | head -n -1 > {params.output_prefix}_unrel_ids.txt
        plink2 --pfile {params.input_prefix} --keep {params.output_prefix}_unrel_ids.txt --make-pgen --out {params.output_prefix}

    elif [[ "{params.method}" == "primus" || "{params.method}" == "2" ]]; then
        echo "PRIMUS ESTIMATION"
        mkdir -p {params.output_prefix}_primus_tmp
        bash {params.scripts_dir}/run_primus.sh {params.input_prefix} {params.output_prefix}_primus_tmp {params.ref_path}
        plink2 --bfile {params.output_prefix}_primus_tmp/unrelated --make-grm-bin --out {params.output_prefix}_grm
        mv {params.output_prefix}_grm.grm.bin {output.grm}
        mv {params.output_prefix}_grm.grm.id {output.grmid}
        mv {params.output_prefix}_grm.grm.N.bin {output.grmN}
        plink2 --bfile {params.output_prefix}_primus_tmp/unrelated --make-pgen --out {params.output_prefix}
        rm -rf {params.output_prefix}_primus_tmp

    else
        echo "ASSUMING UNRELATED: no relatedness method specified"
        cp {input.pgen} {output.pgen}
        cp {input.pvar} {output.pvar}
        cp {input.psam} {output.psam}
        plink2 --pfile {params.input_prefix} --make-grm-bin --out {params.output_prefix}_grm
        mv {params.output_prefix}_grm.grm.bin {output.grm}
        mv {params.output_prefix}_grm.grm.id {output.grmid}
        mv {params.output_prefix}_grm.grm.N.bin {output.grmN}
    fi

    """
