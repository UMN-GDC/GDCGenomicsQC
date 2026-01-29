rule checkRelatedness:
    resources:
        # nodes=1 is usually the default, but can be specified if needed
        nodes = 1,
        # mem=32GB translated to MB
        mem_mb = 128000,
        runtime =720 
    input:
        bed = lambda wildcards: get_input_by_stage(wildcards) + ".LDpruned.bed",
        bim = lambda wildcards: get_input_by_stage(wildcards) + ".LDpruned.bim",
        fam = lambda wildcards: get_input_by_stage(wildcards) + ".LDpruned.fam"
    output:
        bed =   os.path.join(config['OUT_DIR'], "{stage}/unrelated.bed"),
        bim =   os.path.join(config['OUT_DIR'], "{stage}/unrelated.bim"),
        fam =   os.path.join(config['OUT_DIR'], "{stage}/unrelated.fam"),
    params:
        output_dir = os.path.join(config['OUT_DIR'], "{stage}"),
        input_prefix = lambda wildcards, input: input.bed[:-4],
        method = config['relatedness']['method'],
        datatype= "full",
        ref= config["REF"]
    shell: """
    
    echo "Estimating genetic relatedness"
    echo {params.method}

    if [[ {params.method} == king || {params.method} == 1 ]]; then
        echo "KING ESTIMATION"
        bash scripts/run_king.sh {params.output_dir} {params.input_prefix} 1
    elif [[ "{params.method}" == "primus" || "{params.method}" == "2" ]]; then
        echo "PRIMUS ESTIMATION"
        bash scripts/run_primus.sh {params.input_prefix} {params.output_dir} {params.ref} {params.datatype}
    # elif [[ "{params.method}" == "0"]]; then
    else
        echo "ASSUMING UNRELATED SAMPLE, since no method of estimating relatedness specified"
        cp {params.input_prefix}.bed {output.bed}
        cp {params.input_prefix}.bim {output.bim}
        cp {params.input_prefix}.fam {output.fam}
    fi


    """

