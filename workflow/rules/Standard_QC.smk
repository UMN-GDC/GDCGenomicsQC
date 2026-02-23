rule Standard_QC:
    container: "oras://ghcr.io/coffm049/gdcgnomicsqc/plink:latest"
    conda: "../../envs/ancNreport.yml"
    threads: 8
    resources:
        # nodes=1 is usually the default, but can be specified if needed
        nodes = 1,
        # mem=32GB translated to MB
        mem_mb = 32000,
        runtime =60,
    input:
        bed = os.path.join(config['OUT_DIR'], "{stage}/initialFilter.bed"),
        bim = os.path.join(config['OUT_DIR'], "{stage}/initialFilter.bim"),
        fam = os.path.join(config['OUT_DIR'], "{stage}/initialFilter.fam"),
        LDbed = os.path.join(config['OUT_DIR'], "{stage}/initialFilter.LDpruned.bed"),
        LDbim = os.path.join(config['OUT_DIR'], "{stage}/initialFilter.LDpruned.bim"),
        LDfam = os.path.join(config['OUT_DIR'], "{stage}/initialFilter.LDpruned.fam"),
    output:
        bed =   os.path.join(config['OUT_DIR'], "{stage}/standardFiltered.bed"),
        bim =   os.path.join(config['OUT_DIR'], "{stage}/standardFiltered.bim"),
        fam =   os.path.join(config['OUT_DIR'], "{stage}/standardFiltered.fam"),
        LDbed = os.path.join(config['OUT_DIR'], "{stage}/standardFiltered.LDpruned.bed"),
        LDbim = os.path.join(config['OUT_DIR'], "{stage}/standardFiltered.LDpruned.bim"),
        LDfam = os.path.join(config['OUT_DIR'], "{stage}/standardFiltered.LDpruned.fam"),
        tempDir  = temp(directory(os.path.join(config['OUT_DIR'], "{stage}/intermediates/standard_filter/")))
    params:
        ref= config["REF"],
        output_dir = os.path.join(config['OUT_DIR'], "{stage}"),
        sex_check = config['SEX_CHECK'],
        input_prefix = lambda wildcards, input: input.LDbed[:-4],
        relatedness = config['relatedness']['method']
    shell: """
        echo "Standard QC: Variants and samples filtering"
        echo "Input: {wildcards.stage}/{params.input_prefix}"
        echo "Output: {params.output_dir}/standardFiltered"

        if {params.sex_check} ; then
          echo "Performing Sex check"
          plink --bfile {params.input_prefix} --check-sex --out {params.input_prefix}
          grep 'PROBLEM' {params.input_prefix}.sexcheck | awk '{{print $1,$2}}' > {params.output_dir}/sex_discrepancy.txt
          plink --bfile {params.input_prefix} --remove {params.output_dir}/sex_discrepancy.txt --make-bed --out {params.input_prefix}/pastSex
        else
          echo "Ignoring Sex check"
          mv {params.input_prefix}.bed {params.output_dir}/pastSex.bed
          mv {params.input_prefix}.bim {params.output_dir}/pastSex.bim
          mv {params.input_prefix}.fam {params.output_dir}/pastSex.fam
        fi
        bash scripts/filterStandard.sh {params.output_dir}/pastSex {params.output_dir} {threads}
        """
