SLURM_LOGS = "--job-name=%x --output=logs/%x_%j.out --error=logs/%x_%j.err"
rule Initial_QC:
    resources:
        # nodes=1 is usually the default, but can be specified if needed
        nodes = 1,
        # mem=32GB translated to MB
        mem_mb = 32000,
        runtime =60 
    input:
        bed = f"{config['INPUT_FILE']}.bed",
        bim = f"{config['INPUT_FILE']}.bim",
        fam = f"{config['INPUT_FILE']}.fam"
    output:
        # List all files that PLINK will actually create
        bed = f"{config['OUT_DIR']}/Initial_QC/final.bed",
        bim = f"{config['OUT_DIR']}/Initial_QC/final.bim",
        fam = f"{config['OUT_DIR']}/Initial_QC/final.fam",
    params :
        input_prefix = lambda wildcards, input: input.bed[:-4],
        output_dir = config['OUT_DIR'],
    shell: """
    echo "Initial QC: Variants and samples filtering"
    echo "Input: {input}"
    echo "Output: {output}"

    bash scripts/initial_QC.sh {params.input_prefix} {params.output_dir}/Initial_QC
    """

