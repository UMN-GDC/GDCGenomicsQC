rule checkRelatedness:
    container: "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
    conda: "../../envs/ancNreport.yml"
    resources:
        # nodes=1 is usually the default, but can be specified if needed
        nodes = 1,
        # mem=32GB translated to MB
        mem_mb = 128000,
        runtime =720,
    input:
        bed = OUT_DIR / "{subset}" / "initialFilter.LDpruned.bed",
        bim = OUT_DIR / "{subset}" / "initialFilter.LDpruned.bim",
        fam = OUT_DIR / "{subset}" / "initialFilter.LDpruned.fam",
    output:
        bed = OUT_DIR / "{subset}" / "unrelated.bed",
        bim = OUT_DIR / "{subset}" / "unrelated.bim",
        fam = OUT_DIR / "{subset}" / "unrelated.fam",
    params:
        output_dir = OUT_DIR / "{subset}",
        input_prefix = lambda wildcards, input: input.bed[:-4],
        method = config['relatedness']['method'],
        king_cutoff = config['relatedness'].get('king_cutoff', 0.0884),
        ref= REF
    shell: """
    
    echo "Estimating genetic relatedness"
    echo {params.method}

    if [[ {params.method} == king || {params.method} == 1 ]]; then
        echo "KING ESTIMATION with cutoff {params.king_cutoff}"
        bash scripts/run_king.sh {params.output_dir} {params.input_prefix} 1 {params.king_cutoff}
    elif [[ "{params.method}" == "primus" || "{params.method}" == "2" ]]; then
        # [ ] NEeds to be fixed
        echo "PRIMUS ESTIMATION"
        bash scripts/run_primus.sh {params.input_prefix} {params.output_dir} {params.ref}
    # elif [[ "{params.method}" == "0"]]; then
    else
        echo "ASSUMING UNRELATED SAMPLE, since no method of estimating relatedness specified"
        cp {input.bed} {output.bed}
        cp {input.bim} {output.bim}
        cp {input.fam} {output.fam}
    fi


    """

