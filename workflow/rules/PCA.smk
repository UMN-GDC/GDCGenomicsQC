SLURM_LOGS = "--job-name=%x --output=logs/%x_%j.out --error=logs/%x_%j.err"
rule PCA:
    container: "images/my_tool.sif"  # Works for .sif, .img, or docker://
    resources:
        # nodes=1 is usually the default, but can be specified if needed
        nodes = 1,
        # mem=32GB translated to MB
        mem_mb = 32000,
        runtime = 2880
    input:
        bed = f"{config['OUT_DIR']}/relatedness/unrelated.bed",
        bim = f"{config['OUT_DIR']}/relatedness/unrelated.bim",
        fam = f"{config['OUT_DIR']}/relatedness/unrelated.fam",
    output:
        # List all files that PLINK will actually create
        eigen = f"{config['OUT_DIR']}/PCA/merged_dataset_pca.eigenvec",
    params:
        method = config['relatedness']["method"],
        out_dir = f"{config['OUT_DIR']}/PCA",
        input_prefix = f"{config['OUT_DIR']}/relatedness/unrelated",
        ref= config["REF"]
    shell: """
    echo "PCA: "

    echo {params.method}

    if [[ {params.method} == king || {params.method} == 1 || {params.method} == primus || {params.method} == 2 ]]; then
        echo "PCAIR"
    else
        echo "Standard PCA since no method of relatedness estimation included"
        bash scripts/run_pca.sh {params.input_prefix} {params.out_dir} {params.ref}
    fi
    """
