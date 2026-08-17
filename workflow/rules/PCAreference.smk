rule runPcaOnReferencePanel:
    log:
        OUT_DIR / "logs" / "runPcaOnReferencePanel.log",
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
    conda:
        "../../envs/ancNreport.yml"
    threads: 8
    resources:
        nodes=1,
        mem_mb=32000,
        runtime=2880,
    input:
        pgen=OUT_DIR / "full" / "initialFilter.pgen",
        pvar=OUT_DIR / "full" / "initialFilter.pvar",
        psam=OUT_DIR / "full" / "initialFilter.psam",
        ref_bed=lambda wildcards: config.get("ancestry", {}).get(
            "reference_panel_prefix",
            str(REF / "1000G_GRCh38" / "1000G.ensembl.105.with.rsid.gender"),
        )
        + ".bed",
        ref_bim=lambda wildcards: config.get("ancestry", {}).get(
            "reference_panel_prefix",
            str(REF / "1000G_GRCh38" / "1000G.ensembl.105.with.rsid.gender"),
        )
        + ".bim",
        ref_fam=lambda wildcards: config.get("ancestry", {}).get(
            "reference_panel_prefix",
            str(REF / "1000G_GRCh38" / "1000G.ensembl.105.with.rsid.gender"),
        )
        + ".fam",
    output:
        eigen=OUT_DIR / "01-globalAncestry" / "ref.eigenvec.allele",
        projected=OUT_DIR / "01-globalAncestry" / "sampleRefPCscores.sscore",
        projectedref=OUT_DIR / "01-globalAncestry" / "refRefPCscores.sscore",
        tempDir=temp(directory(OUT_DIR / "01-globalAncestry" / "intermediates")),
    params:
        input_prefix=OUT_DIR / "full" / "initialFilter",
        dir=str(OUT_DIR / "01-globalAncestry"),
        ref_prefix=config.get(
            "ancestry",
            {},
        ).get(
            "reference_panel_prefix",
            str(REF / "1000G_GRCh38" / "1000G.ensembl.105.with.rsid.gender"),
        ),
        overlap_snp_list=config.get(
            "ancestry",
            {},
        ).get(
            "overlap_snp_list",
            str(
                OUT_DIR
                / "01-globalAncestry"
                / "intermediates"
                / "abcd_1000ggrch38_overlap.snplist"
            ),
        ),
        overlap_from_bim=config.get("ancestry", {}).get("overlap_from_bim", ""),
        pca_npc=int(config.get("ancestry", {}).get("pca_npc", 10)),
        score_cols=lambda wildcards: f"6-{5 + int(config.get('ancestry', {}).get('pca_npc', 10))}",
    shell:
        r"""
        set -euo pipefail

        echo "PCA:"
        mkdir -p {output.tempDir}

        if [ ! -s "{params.overlap_snp_list}" ]; then
            if [ -z "{params.overlap_from_bim}" ]; then
                echo "No ancestry.overlap_from_bim configured and overlap SNP list missing."
                exit 1
            fi

            awk '{{print $2}}' "{params.overlap_from_bim}" | sort -u > {params.dir}/intermediates/study_rsid.tmp
            awk '{{print $2}}' "{params.ref_prefix}.bim" | sort -u > {params.dir}/intermediates/ref_rsid.tmp
            comm -12 {params.dir}/intermediates/study_rsid.tmp {params.dir}/intermediates/ref_rsid.tmp > "{params.overlap_snp_list}"
            rm -f {params.dir}/intermediates/study_rsid.tmp {params.dir}/intermediates/ref_rsid.tmp
        fi

        if [ ! -s "{params.overlap_snp_list}" ]; then
            echo "Overlap SNP list is empty: {params.overlap_snp_list}"
            exit 1
        fi

        plink2 --bfile {params.ref_prefix} \
               --freq counts \
               --threads {threads} \
               --extract {params.overlap_snp_list} \
               --pca approx {params.pca_npc} allele-wts vcols=chrom,ref,alt \
               --out {params.dir}/ref \
               --allow-no-sex

        echo "Project sample onto the reference PCs."
        plink2 --pfile {params.input_prefix} \
               --read-freq {params.dir}/ref.acount \
               --score {params.dir}/ref.eigenvec.allele 2 5 header-read \
               --score-col-nums {params.score_cols} \
               --extract {params.overlap_snp_list} \
               --out {params.dir}/sampleRefPCscores

        echo "Project ref onto the reference PCs."
        plink2 --bfile {params.ref_prefix} \
               --read-freq {params.dir}/ref.acount \
               --score {params.dir}/ref.eigenvec.allele 2 5 header-read \
               --score-col-nums {params.score_cols} \
               --extract {params.overlap_snp_list} \
               --out {params.dir}/refRefPCscores
        """
