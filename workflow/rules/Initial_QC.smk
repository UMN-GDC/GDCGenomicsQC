rule Initial_QC:
    container: "images/baseImage.sif"
    conda: "../../envs/qcEnvironment.yml"
    resources:
        nodes = 1,
        mem_mb = 32000,
        runtime = 60
    input:
        bed = lambda wildcards: get_input_by_stage(wildcards) + ".bed",
        bim = lambda wildcards: get_input_by_stage(wildcards) + ".bim",
        fam = lambda wildcards: get_input_by_stage(wildcards) + ".fam"
    output:
        bed =   os.path.join(config['OUT_DIR'], "{stage}/initialFilter.bed"),
        bim =   os.path.join(config['OUT_DIR'], "{stage}/initialFilter.bim"),
        fam =   os.path.join(config['OUT_DIR'], "{stage}/initialFilter.fam"),
        LDbed = os.path.join(config['OUT_DIR'], "{stage}/initialFilter.LDpruned.bed"),
        LDbim = os.path.join(config['OUT_DIR'], "{stage}/initialFilter.LDpruned.bim"),
        LDfam = os.path.join(config['OUT_DIR'], "{stage}/initialFilter.LDpruned.fam"),
        tempDir  = temp(directory(os.path.join(config['OUT_DIR'], "{stage}/intermediates/")))
    params:
        # input.bed is a list of one file, we take index 0
        input_prefix = lambda wildcards, input: input.bed[:-4],
        output_prefix = lambda wildcards: os.path.join(config['OUT_DIR'], wildcards.stage)
    shell: """
    echo "Processing Stage: {wildcards.stage}"
    echo "Source Prefix: {params.input_prefix}"
    echo "Output Prefix: {params.output_prefix}/initialFilter"

    bash scripts/initialFilter.sh {params.input_prefix} {params.output_prefix}
    """
