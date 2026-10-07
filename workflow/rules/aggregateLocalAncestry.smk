# Shared aggregation of per-chromosome local ancestry into global proportions.
# Consumes chr{CHR}.lai.msp.tsv + chr{CHR}.lai.fb.tsv produced by EITHER the
# RFMIX estimator (rules/RFMIX.smk) or the Gnomix estimator (rules/gnomix.smk),
# both of which normalize to the same .lai file conventions.
rule aggregateLocalAncestryResults:
    log:
        OUT_DIR / "logs" / "aggregateLocalAncestryResults.log",
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
    conda:
        "../../envs/ancNreport.yml"
    envmodules: *([config.get("R_module")] if config.get("R_module") else [])
    threads: 4
    resources:
        nodes=1,
        mem_mb=32000,
        runtime=60,
    input:
        msp=expand(
            OUT_DIR / "02-localAncestry" / "chr{CHR}.lai.msp.tsv", CHR=LOCAL_ANCESTRY_CHROMOSOMES
        ),
        fb=expand(OUT_DIR / "02-localAncestry" / "chr{CHR}.lai.fb.tsv", CHR=LOCAL_ANCESTRY_CHROMOSOMES),
    output:
        mat=OUT_DIR / "02-localAncestry" / "ancestry_full.txt",
    params:
        script=workflow.source_path("../scripts/aggregateLocalAncestry.R"),
        out_dir=OUT_DIR,
    shell: """
        Rscript {params.script} {params.out_dir}
"""
