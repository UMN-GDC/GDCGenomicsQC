rule checkRelatednessExtractUnrelated:
    log:
        OUT_DIR / "logs" / "checkRelatednessExtractUnrelated_{subset}.log",
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
    conda:
        "../../envs/ancNreport.yml"
    envmodules: lambda wildcards: [config["plink_module"]] if config.get("plink_module") else []
    resources:
        nodes=1,
        mem_mb=128000,
        runtime=720,
    threads: 16
    input:
        pgen=OUT_DIR / "{subset}" / "f1.b38.ldpruned.pgen",
        pvar=OUT_DIR / "{subset}" / "f1.b38.ldpruned.pvar",
        psam=OUT_DIR / "{subset}" / "f1.b38.ldpruned.psam",
    output:
        pgen=OUT_DIR / "{subset}" / "f1.b38.ldpruned.unrelated.pgen",
        pvar=OUT_DIR / "{subset}" / "f1.b38.ldpruned.unrelated.pvar",
        psam=OUT_DIR / "{subset}" / "f1.b38.ldpruned.unrelated.psam",
        unrels=OUT_DIR / "{subset}" / "f1.b38.ldpruned.unrelated_ids.txt",
        grm=OUT_DIR / "{subset}" / "f1.b38.ldpruned.unrelated.grm.bin",
        grmid=OUT_DIR / "{subset}" / "f1.b38.ldpruned.unrelated.grm.id",
        grmN=OUT_DIR / "{subset}" / "f1.b38.ldpruned.unrelated.grm.N.bin",
    params:
        king_cutoff=config.get("relatedness", {}).get("king_cutoff", 0.0884),
        method=config.get("relatedness", {}).get("method", "king"),
        scripts_dir=SCRIPTS_DIR,
        input_prefix=lambda wildcards, input: input.pgen[:-5],
        output_prefix=lambda wildcards, output: output.pgen[:-5],
        ref_path=config.get("REF", "/path/to/ref"),
    shell:
        """

    echo "Estimating genetic relatedness"
    echo "Method: {params.method}"

    if [[ "{params.method}" == "king" || "{params.method}" == "1" ]]; then
        echo "KING ESTIMATION using PLINK2 with cutoff {params.king_cutoff}"
        plink2 --pfile {params.input_prefix} \
            --make-grm-bin \
            --threads {threads} \
            --king-cutoff {params.king_cutoff} \
            --make-king \
            --out {params.output_prefix}_grm
        mv {params.output_prefix}_grm.grm.bin {output.grm}
        mv {params.output_prefix}_grm.grm.id {output.grmid}
        mv {params.output_prefix}_grm.grm.N.bin {output.grmN}
        awk -v cutoff={params.king_cutoff} '$5 < cutoff {{print $1, $2}}' {params.output_prefix}_grm.king.cutoff.in.id | head -n -1 > {output.unrels}
        plink2 --pfile {params.input_prefix} --keep {output.unrels} --make-pgen --out {params.output_prefix}

    elif [[ "{params.method}" == "primus" || "{params.method}" == "2" ]]; then
        echo "PRIMUS ESTIMATION"
        mkdir -p {params.output_prefix}_primus_tmp
        bash {params.scripts_dir}/run_primus.sh {params.input_prefix} {params.output_prefix}_primus_tmp {params.ref_path}
        plink2 --bfile {params.output_prefix}_primus_tmp/unrelated --make-grm-bin --out {params.output_prefix}_grm --threads {threads}
        mv {params.output_prefix}_grm.grm.bin {output.grm}
        mv {params.output_prefix}_grm.grm.id {output.grmid}
        mv {params.output_prefix}_grm.grm.N.bin {output.grmN}
        plink2 --bfile {params.output_prefix}_primus_tmp/unrelated --make-pgen --out {params.output_prefix} --threads {threads}
        rm -rf {params.output_prefix}_primus_tmp

    else
        echo "ASSUMING UNRELATED: no relatedness method specified"
        cp {input.pgen} {output.pgen}
        cp {input.pvar} {output.pvar}
        cp {input.psam} {output.psam}
        plink2 --pfile {params.input_prefix} --make-grm-bin --out {params.output_prefix}_grm --threads {threads}
        mv {params.output_prefix}_grm.grm.bin {output.grm}
        mv {params.output_prefix}_grm.grm.id {output.grmid}
        mv {params.output_prefix}_grm.grm.N.bin {output.grmN}
    fi

    """


rule ldPruneUnrelated:
    log:
        OUT_DIR / "logs" / "ldPruneUnrelated_{subset}.log",
    container:
        "docker://gfanz/plink2:latest"
    conda:
        "../../envs/ancNreport.yml"
    envmodules: lambda wildcards: [config["plink_module"]] if config.get("plink_module") else []
    threads: 8
    resources:
        nodes=1,
        mem_mb=32000,
        runtime=60,
    input:
        pgen=OUT_DIR / "{subset}" / "f1.b38.ldpruned.unrelated.pgen",
        pvar=OUT_DIR / "{subset}" / "f1.b38.ldpruned.unrelated.pvar",
        psam=OUT_DIR / "{subset}" / "f1.b38.ldpruned.unrelated.psam",
    output:
        pgen=OUT_DIR / "{subset}" / "f1.b38.ldpruned.unrelated.ldpruned.pgen",
        pvar=OUT_DIR / "{subset}" / "f1.b38.ldpruned.unrelated.ldpruned.pvar",
        psam=OUT_DIR / "{subset}" / "f1.b38.ldpruned.unrelated.ldpruned.psam",
    params:
        input_prefix=lambda wildcards, input: input.pgen[:-5],
    shell:
        """
        echo "LD pruning unrelated samples"
        plink2 --pfile {params.input_prefix} --indep-pairwise 50 5 0.2 --out {params.input_prefix}_indep --threads {threads}
        plink2 --pfile {params.input_prefix} --extract {params.input_prefix}_indep.prune.in --make-pgen --out {params.input_prefix}.ldpruned --threads {threads}
        """
