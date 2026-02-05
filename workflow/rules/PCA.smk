rule PCA:
    container: "../envs/plink.sif"
    conda: "../../envs/plink.yml"
    threads: 8
    resources:
        nodes = 1,
        mem_mb = 32000,
        runtime = 2880,
    input:
        bed = f"{config['OUT_DIR']}/02-relatedness/standardFiltered.LDpruned.bed",
        bim = f"{config['OUT_DIR']}/02-relatedness/standardFiltered.LDpruned.bim",
        fam = f"{config['OUT_DIR']}/02-relatedness/standardFiltered.LDpruned.fam",
    output:
        # List all files that PLINK will actually create
        eigen = f"{config['OUT_DIR']}/04-globalAncestry/merged_dataset_pca.eigenvec",
        tempDir  = temp(directory(f"{config['OUT_DIR']}/04-globalAncestry/intermediates/"))
    params:
        method = config['relatedness']["method"],
        grm = config['relatedness']["method"],
        out_dir = f"{config['OUT_DIR']}/04-globalAncestry",
        input_prefix = f"{config['OUT_DIR']}/02-relatedness/standardFiltered.LDpruned",
        ref= config["REF"]
    shell: """
    echo "PCA: "

    echo {params.method}

    if [[ {params.method} == king || {params.method} == 1 || {params.method} == primus || {params.method} == 2 ]]; then
        echo "PCAIR"
    else
        echo "Standard PCA since no method of relatedness estimation included"
        bash scripts/run_pca.sh {params.input_prefix} {params.out_dir} {params.ref} {params.grm}
    fi
    """
