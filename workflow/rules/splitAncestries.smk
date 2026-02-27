rule splitAncestries:
    conda: "../../envs/ancNreport.yml"
    resources:
        nodes = 1,
        mem_mb = 8000,
        runtime = 15,
    input:
        bed = OUT_DIR / "00-raw" / "data.bed",
        bim = OUT_DIR / "00-raw" / "data.bim",
        fam = OUT_DIR / "00-raw" / "data.fam",
        anc = OUT_DIR / "04-globalAncestry" / "latentDistantRelatedness.tsv",
    output:
        bed =  OUT_DIR / "{ANCESTRY}" / "data.bed",
        bim = OUT_DIR / "{ANCESTRY}" / "data.bim",
        fam = OUT_DIR / "{ANCESTRY}" / "data.fam",
    params:
        outdir = config['OUT_DIR'],
        method = "pc"

    shell: """
    echo "Splitting ancestries"

    plink2 \
      --bed {input.bed} \
      --bim {input.bim} \
      --fam {input.fam} \
      --covar {input.anc} \
      --keep-if 'pc_label == {wildcards.ANCESTRY}' \
      --out {params.outdir}/{wildcards.ANCESTRY}/data \
      --make-bed
    """


