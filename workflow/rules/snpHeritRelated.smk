# SNP heritability estimation on samples including close relatives (fastGWA
# and/or PredLMM-ACE). Can also be run on an unrelated homogeneous subset
# (e.g. AFR) by setting snpHeritRelated.related: false and subset: "AFR".

SNP_HERIT_RELATED_CONFIG = config.get("snpHeritRelated", {})
SNP_HERIT_RELATED_ACTIVE = bool(SNP_HERIT_RELATED_CONFIG.get("pheno"))
SNP_HERIT_RELATED_METHODS = SNP_HERIT_RELATED_CONFIG.get("methods", ["fastgwa", "predlmm_ace"])
SNP_HERIT_RELATED_METHODS = [str(m).lower() for m in SNP_HERIT_RELATED_METHODS]
SNP_HERIT_RELATED_OUTDIR = SNP_HERIT_RELATED_CONFIG.get("output_dir", "03-snpHeritability")
SNP_HERIT_RELATED_SUBSET = SNP_HERIT_RELATED_CONFIG.get("subset")

# Image built from envs/predlmmAce.def (see envs/README.md). Quoted so the
# oras:// prefix is not mistaken for a path; used together with the conda
# directive so the cluster can still build envs/predlmmAce.yml directly.
PREDLMM_ACE_CONTAINER = "oras://ghcr.io/coffm049/gdcgenomicsqc/predlmmace:latest"

if SNP_HERIT_RELATED_CONFIG:
    valid_methods = {"fastgwa", "predlmm_ace"}
    bad = set(SNP_HERIT_RELATED_METHODS) - valid_methods
    if bad:
        raise ValueError(f"snpHeritRelated.methods must be in {valid_methods}, got {bad}")
    if SNP_HERIT_RELATED_ACTIVE:
        has_out = bool(SNP_HERIT_RELATED_CONFIG.get("out"))
        has_grm = bool(SNP_HERIT_RELATED_CONFIG.get("grm_prefix"))
        has_pca = bool(SNP_HERIT_RELATED_CONFIG.get("pca_input"))
        has_bed = bool(SNP_HERIT_RELATED_CONFIG.get("bed_prefix"))
        if has_out or has_grm or has_pca or has_bed:
            if not (has_out and has_grm and has_pca):
                raise ValueError(
                    "snpHeritRelated: external mode requires all of (out, grm_prefix, pca_input)"
                    " and bed_prefix when fastgwa is requested"
                )
            if "fastgwa" in SNP_HERIT_RELATED_METHODS and not has_bed:
                raise ValueError(
                    "snpHeritRelated: fastgwa requires bed_prefix in external mode"
                )
            SNP_HERIT_RELATED_EXTERNAL = True
        else:
            SNP_HERIT_RELATED_EXTERNAL = False


def _shr_dense_prefix(subset=None):
    """Dense-GRM path prefix: related=True -> full sample set (incl. relatives);
    related=False -> unrelated LD-pruned set within the subset."""
    if SNP_HERIT_RELATED_CONFIG.get("related", True) is False:
        name = "f1.b38.ldpruned.unrelated.ldpruned"
    else:
        name = "f1.b38.ldpruned"
    return OUT_DIR / subset / name if subset else None


# Output directives may not be lambdas ("Only input files can be specified as
# functions"), so the dense-set base name is resolved once at parse time and
# the {subset} wildcard is carried by the plain path strings below.
_SHR_BASE = (
    "f1.b38.ldpruned.unrelated.ldpruned"
    if SNP_HERIT_RELATED_CONFIG.get("related", True) is False
    else "f1.b38.ldpruned"
)


if SNP_HERIT_RELATED_ACTIVE:

    rule convertToBedForRelatedHerit:
        log:
            OUT_DIR / "logs" / "convertToBedForRelatedHerit_{subset}.log",
        container:
            "docker://gfanz/plink2:latest",
        conda:
            "../../envs/ancNreport.yml",
        envmodules: *([config.get("plink_module")] if config.get("plink_module") else []),
        threads: 8,
        resources:
            nodes=1,
            mem_mb=32000,
            runtime=60,
        input:
            pgen=lambda w: OUT_DIR / w.subset / (_SHR_BASE + ".pgen"),
            pvar=lambda w: OUT_DIR / w.subset / (_SHR_BASE + ".pvar"),
            psam=lambda w: OUT_DIR / w.subset / (_SHR_BASE + ".psam"),
        output:
            bed=OUT_DIR / "{subset}" / (_SHR_BASE + ".bed"),
            bim=OUT_DIR / "{subset}" / (_SHR_BASE + ".bim"),
            fam=OUT_DIR / "{subset}" / (_SHR_BASE + ".fam"),
        params:
            in_prefix=lambda w, input: str(input.pgen)[:-5],
            out_prefix=lambda w, output: str(output.bed)[:-4],
        shell:
            """
            plink2 --pfile {params.in_prefix} --make-bed --out {params.out_prefix} --threads {threads}
            """

    rule makeSparseGrmForRelatedHerit:
        log:
            OUT_DIR / "logs" / "makeSparseGrmForRelatedHerit_{subset}.log",
        conda:
            "../../envs/predlmmAce.yml",
        container:
            PREDLMM_ACE_CONTAINER,
        threads: 8,
        resources:
            nodes=1,
            mem_mb=32000,
            runtime=360,
        input:
            grm_bin=lambda w: OUT_DIR / w.subset / (_SHR_BASE + ".grm.bin"),
        output:
            sparse=OUT_DIR / "{subset}" / (_SHR_BASE + "_sparsegrm.grm.sp"),
        params:
            grm_prefix=lambda w: str(_shr_dense_prefix(w.subset)),
            out_prefix=lambda w, output: str(output.sparse)[: -len(".grm.sp")],
            cutoff=SNP_HERIT_RELATED_CONFIG.get("sparse_cutoff", 0.05),
        shell:
            """
            gcta --grm {params.grm_prefix} --make-bK-sparse {params.cutoff} --out {params.out_prefix} --thread-num {threads}
            """

    rule estimateRelatedHeritFastGWA:
        log:
            OUT_DIR / "logs" / "estimateRelatedHeritFastGWA_{subset}.log",
        conda:
            "../../envs/predlmmAce.yml",
        container:
            PREDLMM_ACE_CONTAINER,
        threads: 16,
        resources:
            nodes=1,
            mem_mb=64000,
            runtime=1440,
        input:
            bed=lambda w: OUT_DIR / w.subset / (_SHR_BASE + ".bed"),
            bim=lambda w: OUT_DIR / w.subset / (_SHR_BASE + ".bim"),
            fam=lambda w: OUT_DIR / w.subset / (_SHR_BASE + ".fam"),
            sparse=lambda w: OUT_DIR / w.subset / (_SHR_BASE + "_sparsegrm.grm.sp"),
            eigenvec=OUT_DIR / "{subset}" / "internal_pca_plink2.eigenvec",
        output:
            fastgwa=OUT_DIR / "{subset}" / SNP_HERIT_RELATED_OUTDIR / "fastgwa.fastGWA",
            fastgwa_log=OUT_DIR / "{subset}" / SNP_HERIT_RELATED_OUTDIR / "fastgwa.log",
        params:
            pheno=SNP_HERIT_RELATED_CONFIG["pheno"],
            covar=SNP_HERIT_RELATED_CONFIG.get("covar"),
            mpheno=SNP_HERIT_RELATED_CONFIG.get("mpheno", 1),
            bed_prefix=lambda w: str(_shr_dense_prefix(w.subset)),
            sparse_prefix=lambda w: str(_shr_dense_prefix(w.subset)) + "_sparsegrm",
            gcta_prefix=lambda w: str(OUT_DIR / w.subset / SNP_HERIT_RELATED_OUTDIR / "fastgwa"),
            covar_flag="--covar " + str(SNP_HERIT_RELATED_CONFIG["covar"]) if SNP_HERIT_RELATED_CONFIG.get("covar") else "",
            scripts_dir=SCRIPTS_DIR,
        shell:
            """
            mkdir -p "$(dirname {output.fastgwa})"
            gcta --bfile {params.bed_prefix} \\
                --grm-sparse {params.sparse_prefix} \\
                --fastGWA-mlm \\
                --pheno {params.pheno} \\
                --mpheno {params.mpheno} \\
                --qcovar {input.eigenvec} \\
                {params.covar_flag} \\
                --out {params.gcta_prefix}
            """

    rule finalizeRelatedHeritFastGWA:
        log:
            OUT_DIR / "logs" / "finalizeRelatedHeritFastGWA_{subset}.log",
        conda:
            "../../envs/predlmmAce.yml",
        container:
            PREDLMM_ACE_CONTAINER,
        threads: 1,
        resources:
            nodes=1,
            mem_mb=8000,
            runtime=30,
        input:
            fastgwa=OUT_DIR / "{subset}" / SNP_HERIT_RELATED_OUTDIR / "fastgwa.fastGWA",
            fastgwa_log=OUT_DIR / "{subset}" / SNP_HERIT_RELATED_OUTDIR / "fastgwa.log",
        output:
            estimates=OUT_DIR / "{subset}" / SNP_HERIT_RELATED_OUTDIR / "fastgwa_h2.csv",
        params:
            scripts_dir=SCRIPTS_DIR,
        shell:
            """
            python {params.scripts_dir}/parse_fastgwa.py --log {input.fastgwa_log} --out {output.estimates}
            """

    rule estimateRelatedHeritPredLMMAce:
        log:
            OUT_DIR / "logs" / "estimateRelatedHeritPredLMMAce_{subset}.log",
        conda:
            "../../envs/predlmmAce.yml",
        container:
            PREDLMM_ACE_CONTAINER,
        threads: 8,
        resources:
            nodes=1,
            mem_mb=32000,
            runtime=1440,
        input:
            grm_bin=lambda w: OUT_DIR / w.subset / (_SHR_BASE + ".grm.bin"),
            grm_id=lambda w: OUT_DIR / w.subset / (_SHR_BASE + ".grm.id"),
            grm_Nbin=lambda w: OUT_DIR / w.subset / (_SHR_BASE + ".grm.N.bin"),
            eigenvec=OUT_DIR / "{subset}" / "internal_pca_plink2.eigenvec",
        output:
            estimates=OUT_DIR / "{subset}" / SNP_HERIT_RELATED_OUTDIR / "predlmm_ace_h2.csv",
        params:
            pheno=SNP_HERIT_RELATED_CONFIG["pheno"],
            covar=SNP_HERIT_RELATED_CONFIG.get("covar"),
            npc=SNP_HERIT_RELATED_CONFIG.get("npc", 10),
            model=SNP_HERIT_RELATED_CONFIG.get("model", "auto"),
            rank=SNP_HERIT_RELATED_CONFIG.get("rank", 5000),
            family_column=SNP_HERIT_RELATED_CONFIG.get("family_column", "FID"),
            grm_prefix=lambda w: str(_shr_dense_prefix(w.subset)),
            covar_flag="--covar " + str(SNP_HERIT_RELATED_CONFIG["covar"]) if SNP_HERIT_RELATED_CONFIG.get("covar") else "",
            scripts_dir=SCRIPTS_DIR,
        shell:
            """
            mkdir -p "$(dirname {output.estimates})"
            python {params.scripts_dir}/run_predlmm_ace.py \\
                --grm-prefix {params.grm_prefix} \\
                --pheno {params.pheno} \\
                {params.covar_flag} \\
                --eigenvec {input.eigenvec} \\
                --npc {params.npc} \\
                --family-column {params.family_column} \\
                --model {params.model} \\
                --rank {params.rank} \\
                --out {output.estimates}
            """

    if not SNP_HERIT_RELATED_EXTERNAL:

        # DAG mode aggregates across subsets like the other run_ rules
        pass

    # External mode: direct inputs, no {subset} wildcard
    if SNP_HERIT_RELATED_EXTERNAL:

        rule convertToBed_ext:
            conda:
                "../../envs/ancNreport.yml",
            container:
                "docker://gfanz/plink2:latest",
            threads: 8,
            resources:
                nodes=1, mem_mb=32000, runtime=60,
            input:
                bed=SNP_HERIT_RELATED_CONFIG["bed_prefix"] + ".bed",
                bim=SNP_HERIT_RELATED_CONFIG["bed_prefix"] + ".bim",
                fam=SNP_HERIT_RELATED_CONFIG["bed_prefix"] + ".fam",
            output:
                done=SNP_HERIT_RELATED_CONFIG["out"] + "_bed.done",
            shell:
                'touch {output.done}'

        rule makeSparseGrm_ext:
            conda:
                "../../envs/predlmmAce.yml",
            container:
                PREDLMM_ACE_CONTAINER,
            threads: 8,
            resources:
                nodes=1, mem_mb=32000, runtime=360,
            input:
                grm_bin=SNP_HERIT_RELATED_CONFIG["grm_prefix"] + ".grm.bin",
            output:
                sparse=SNP_HERIT_RELATED_CONFIG["grm_prefix"] + "_sparsegrm.grm.sp",
            params:
                grm_prefix=SNP_HERIT_RELATED_CONFIG["grm_prefix"],
                out_prefix=SNP_HERIT_RELATED_CONFIG["grm_prefix"] + "_sparsegrm",
                cutoff=SNP_HERIT_RELATED_CONFIG.get("sparse_cutoff", 0.05),
            shell:
                """
                gcta --grm {params.grm_prefix} --make-bK-sparse {params.cutoff} --out {params.out_prefix} --thread-num {threads}
                """

        rule estimateRelatedHeritFastGWA_ext:
            conda:
                "../../envs/predlmmAce.yml",
            container:
                PREDLMM_ACE_CONTAINER,
            threads: 16,
            resources:
                nodes=1, mem_mb=64000, runtime=1440,
            input:
                bed=SNP_HERIT_RELATED_CONFIG["bed_prefix"] + ".bed",
                bim=SNP_HERIT_RELATED_CONFIG["bed_prefix"] + ".bim",
                fam=SNP_HERIT_RELATED_CONFIG["bed_prefix"] + ".fam",
                sparse=SNP_HERIT_RELATED_CONFIG["grm_prefix"] + "_sparsegrm.grm.sp",
                eigenvec=SNP_HERIT_RELATED_CONFIG["pca_input"],
            output:
                fastgwa=SNP_HERIT_RELATED_CONFIG["out"] + "_fastgwa.fastGWA",
                fastgwa_log=SNP_HERIT_RELATED_CONFIG["out"] + "_fastgwa.log",
            params:
                pheno=SNP_HERIT_RELATED_CONFIG["pheno"],
                covar=SNP_HERIT_RELATED_CONFIG.get("covar"),
                mpheno=SNP_HERIT_RELATED_CONFIG.get("mpheno", 1),
                sparse_prefix=SNP_HERIT_RELATED_CONFIG["grm_prefix"] + "_sparsegrm",
                out_prefix=SNP_HERIT_RELATED_CONFIG["out"] + "_fastgwa",
                bed_prefix=SNP_HERIT_RELATED_CONFIG["bed_prefix"],
                covar_flag="--covar " + str(SNP_HERIT_RELATED_CONFIG["covar"]) if SNP_HERIT_RELATED_CONFIG.get("covar") else "",
            shell:
                """
                gcta --bfile {params.bed_prefix} \\
                    --grm-sparse {params.sparse_prefix} \\
                    --fastGWA-mlm \\
                    --pheno {params.pheno} \\
                    --mpheno {params.mpheno} \\
                    --qcovar {input.eigenvec} \\
                    {params.covar_flag} \\
                    --out {params.out_prefix}
                """

        rule finalizeExtFastGWA:
            conda:
                "../../envs/predlmmAce.yml",
            container:
                PREDLMM_ACE_CONTAINER,
            threads: 1,
            resources:
                nodes=1, mem_mb=8000, runtime=30,
            input:
                fastgwa_log=SNP_HERIT_RELATED_CONFIG["out"] + "_fastgwa.log",
            output:
                estimates=SNP_HERIT_RELATED_CONFIG["out"] + "_fastgwa_h2.csv",
            params:
                scripts_dir=SCRIPTS_DIR,
            shell:
                """
                python {params.scripts_dir}/parse_fastgwa.py --log {input.fastgwa_log} --out {output.estimates}
                """

        rule estimateRelatedHeritPredLMMAce_ext:
            conda:
                "../../envs/predlmmAce.yml",
            container:
                PREDLMM_ACE_CONTAINER,
            threads: 8,
            resources:
                nodes=1, mem_mb=32000, runtime=1440,
            input:
                grm_bin=SNP_HERIT_RELATED_CONFIG["grm_prefix"] + ".grm.bin",
                grm_id=SNP_HERIT_RELATED_CONFIG["grm_prefix"] + ".grm.id",
                grm_Nbin=SNP_HERIT_RELATED_CONFIG["grm_prefix"] + ".grm.N.bin",
                eigenvec=SNP_HERIT_RELATED_CONFIG["pca_input"],
            output:
                estimates=SNP_HERIT_RELATED_CONFIG["out"] + "_predlmm_ace_h2.csv",
            params:
                pheno=SNP_HERIT_RELATED_CONFIG["pheno"],
                covar=SNP_HERIT_RELATED_CONFIG.get("covar"),
                npc=SNP_HERIT_RELATED_CONFIG.get("npc", 10),
                model=SNP_HERIT_RELATED_CONFIG.get("model", "auto"),
                rank=SNP_HERIT_RELATED_CONFIG.get("rank", 5000),
                family_column=SNP_HERIT_RELATED_CONFIG.get("family_column", "FID"),
                grm_prefix=SNP_HERIT_RELATED_CONFIG["grm_prefix"],
                covar_flag="--covar " + str(SNP_HERIT_RELATED_CONFIG["covar"]) if SNP_HERIT_RELATED_CONFIG.get("covar") else "",
                scripts_dir=SCRIPTS_DIR,
            shell:
                """
                python {params.scripts_dir}/run_predlmm_ace.py \\
                    --grm-prefix {params.grm_prefix} \\
                    --pheno {params.pheno} \\
                    {params.covar_flag} \\
                    --eigenvec {input.eigenvec} \\
                    --npc {params.npc} \\
                    --family-column {params.family_column} \\
                    --model {params.model} \\
                    --rank {params.rank} \\
                    --out {output.estimates}
                """
