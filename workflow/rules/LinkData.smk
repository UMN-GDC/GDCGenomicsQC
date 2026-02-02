rule linkData:
    container: "../envs/plink.sif"
    conda: "../../envs/plink.yml"
    input:
        bed = f"{config['INPUT_FILE']}.bed",
        bim = f"{config['INPUT_FILE']}.bim",
        fam = f"{config['INPUT_FILE']}.fam"
    output:
        bed = f"{config['OUT_DIR']}/00-raw/data.bed",
        bim = f"{config['OUT_DIR']}/00-raw/data.bim",
        fam = f"{config['OUT_DIR']}/00-raw/data.fam",
        log = f"{config['OUT_DIR']}/00-raw/data.log",
    params:
        thin = f"{config['thin']}",
        input_prefix = lambda wildcards, input: input.bed[:-4],
        output_prefix = lambda wildcards, output: output.bed[:-4]
    shell: """
        # thinning can be useful for testing out pipeline on large datasets 
        mkdir -p {params.output_prefix}

        if [ "{params.thin}" = "True" ] ;  then
          echo "Thinning data for testing purposes"
          plink2 --bfile {params.input_prefix} \
            --thin-indiv-count 3000 \
            --out {params.output_prefix} \
            --make-bed
        else
          ln -s {input.bed} {output.bed}
          ln -s {input.bim} {output.bim}
          ln -s {input.fam} {output.fam}
        fi
    """
