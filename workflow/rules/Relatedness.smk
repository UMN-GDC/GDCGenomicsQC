rule checkRelatednessExtractUnrelated:
    log:
        OUT_DIR / "logs" / "checkRelatednessExtractUnrelated_{subset}.log",
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
    conda:
        "../../envs/ancNreport.yml"
    threads: 8
    resources:
        nodes=1,
        mem_mb=128000,
        runtime=720,
    input:
        pgen=OUT_DIR / "{subset}" / "standardFilter.pgen",
        pvar=OUT_DIR / "{subset}" / "standardFilter.pvar",
        psam=OUT_DIR / "{subset}" / "standardFilter.psam",
        LDpgen=OUT_DIR / "{subset}" / "standardFilter.LDpruned.pgen",
        LDpvar=OUT_DIR / "{subset}" / "standardFilter.LDpruned.pvar",
        LDpsam=OUT_DIR / "{subset}" / "standardFilter.LDpruned.psam",
    output:
        bed=OUT_DIR / "{subset}" / "unrelated.bed",
        bim=OUT_DIR / "{subset}" / "unrelated.bim",
        fam=OUT_DIR / "{subset}" / "unrelated.fam",
        grm=OUT_DIR / "{subset}" / "unrelated.grm.bin",
        grmid=OUT_DIR / "{subset}" / "unrelated.grm.id",
        grmN=OUT_DIR / "{subset}" / "unrelated.grm.N.bin",
        keep=OUT_DIR / "{subset}" / "unrelated.keep",
        summary=OUT_DIR / "{subset}" / "relatedness_summary.tsv",
        tempDir=temp(
            directory(OUT_DIR / "{subset}" / "intermediates" / "relatedness")
        ),
    params:
        king_cutoff=config.get("relatedness", {}).get("king_cutoff", 0.0884),
        method=config.get("relatedness", {}).get("method", "king"),
        full_prefix=lambda wc, input: str(input.pgen)[:-5],
        ld_prefix=lambda wc, input: str(input.LDpgen)[:-5],
        output_prefix=lambda wc, output: str(output.bed)[:-4],
    shell:
        r"""
        set -euo pipefail

        echo "Estimating genetic relatedness"
        echo "Subset: {wildcards.subset}"
        echo "Method: {params.method}"
        mkdir -p {output.tempDir}

        N_BEFORE=$(( $(wc -l < {input.psam}) - 1 ))

        if [[ "{params.method}" == "king" || "{params.method}" == "1" ]]; then
          echo "KING cutoff: {params.king_cutoff}"

          # Select an unrelated sample set using final QC, LD-pruned variants.
          plink2 --pfile {params.ld_prefix} \
            --king-cutoff {params.king_cutoff} \
            --make-pgen \
            --out {output.tempDir}/king_unrelated \
            --threads {threads}

          awk 'BEGIN {{OFS="\t"}}
            NR==1 {{print "#FID","IID"; next}}
            {{print $1,$2}}' {output.tempDir}/king_unrelated.psam \
            > {output.keep}

        elif [[ "{params.method}" == "0" || "{params.method}" == "none" ]]; then
          echo "Relatedness filtering disabled"
          awk 'BEGIN {{OFS="\t"}}
            NR==1 {{print "#FID","IID"; next}}
            {{print $1,$2}}' {input.psam} > {output.keep}
        else
          echo "ERROR: unsupported relatedness method: {params.method}" >&2
          exit 1
        fi

        N_AFTER=$(( $(wc -l < {output.keep}) - 1 ))
        if [[ "$N_AFTER" -le 0 || "$N_AFTER" -gt "$N_BEFORE" ]]; then
          echo "ERROR: invalid unrelated sample count: $N_AFTER of $N_BEFORE" >&2
          exit 1
        fi

        # Apply the unrelated sample set to all final QC variants.
        plink2 --pfile {params.full_prefix} \
          --keep {output.keep} \
          --make-bed \
          --out {params.output_prefix} \
          --threads {threads}

        # Build the heritability GRM from the full final marker set, not from
        # the LD-pruned relatedness marker set.
        plink2 --bfile {params.output_prefix} \
          --make-grm-bin \
          --out {params.output_prefix} \
          --threads {threads}

        printf "subset\tmethod\tking_cutoff\tn_before\tn_after\tn_removed\n" \
          > {output.summary}
        printf "%s\t%s\t%s\t%d\t%d\t%d\n" \
          "{wildcards.subset}" "{params.method}" "{params.king_cutoff}" \
          "$N_BEFORE" "$N_AFTER" "$((N_BEFORE-N_AFTER))" \
          >> {output.summary}

        echo "Relatedness filtering retained $N_AFTER of $N_BEFORE samples"
        """
