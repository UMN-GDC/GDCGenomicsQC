SLURM_LOGS = "--job-name=%x --output=logs/%x_%j.out --error=logs/%x_%j.err"
rule Relatedness:
    resources:
        # nodes=1 is usually the default, but can be specified if needed
        nodes = 1,
        # mem=32GB translated to MB
        mem_mb = 128000,
        runtime =720 
    input:
        bed = f"{config['OUT_DIR']}/Initial_QC/final.bed",
        fam = f"{config['OUT_DIR']}/Initial_QC/final.fam",
        bim = f"{config['OUT_DIR']}/Initial_QC/final.bim",
    output:
        # List all files that PLINK will actually create
        bed = f"{config['OUT_DIR']}/relatedness/unrelated.bed",
        fam = f"{config['OUT_DIR']}/relatedness/unrelated.fam",
        bim = f"{config['OUT_DIR']}/relatedness/unrelated.bim",
    params:
        output_dir = f"{config['OUT_DIR']}/relatedness",
        input_prefix = f"{config['OUT_DIR']}/Initial_QC/final.LDpruned",
        method = config['relatedness'],
        datatype= "full"
    shell: """
    
    echo "Estimating genetic relatedness"

    if [[ "{params.method}" == "king" || "{params.method}" == "1" ]]; then
        echo "KING ESTIMATION"
        bash scripts/run_king.sh {params.output_dir} {params.input_prefix} 1
    elif [[ "{params.method}" == "primus" || "{params.method}" == "2" ]]; then
        echo "PRIMUS ESTIMATION"
        bash scripts/run_primus.sh {params.input_prefix} {params.output_dir} /projects/standard/gdc/public/Ref {params.datatype}
    # elif [[ {params.method} == 0]]; then
    else
        echo "ASSUMING UNRELATED SAMPLE"
        cp {params.input_prefix}.bed {output.bed}
        cp {params.input_prefix}.bim {output.bim}
        cp {params.input_prefix}.fam {output.fam}
    fi


    """

