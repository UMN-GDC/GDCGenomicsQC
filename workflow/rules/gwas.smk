from pathlib import Path
# Genome-wide association testing on pipeline subsets.
#
# Three methods are supported:
#   standard : PLINK2 --glm (ordinary linear regression with covariates/PCs)
#   fastgwa  : GCTA fastGWA (--fastGWA-lr, or --fastGWA-mlm with a sparse GRM)
#   tractor  : TRACTOR local-ancestry-aware regression using LAI tract output
#
# Subsets follow the same convention as the rest of the pipeline: "full" plus
# every ancestry predicted by the global-ancestry checkpoint (or the explicit
# `gwas.subsets` list). Discovery/target splitting and PRS wiring are described
# in docs/.

GWAS_CONFIG = config.get("gwas", {})
GWAS_ACTIVE = bool(GWAS_CONFIG.get("pheno"))
GWAS_METHODS = [str(m).lower() for m in GWAS_CONFIG.get("methods", ["standard", "fastgwa"])]
GWAS_OUTDIR = GWAS_CONFIG.get("output_dir", "04-gwas")
GWAS_BED_PREFIX = GWAS_CONFIG.get("bed_prefix", "f1.b38")
GWAS_GRM_PREFIX = GWAS_CONFIG.get("grm_prefix", "f1.b38.ldpruned")
GWAS_FASTGWA_MODE = str(GWAS_CONFIG.get("fastgwa_mode", "lr")).lower()
GWAS_PHENO = GWAS_CONFIG.get("pheno")
GWAS_MPHENO = GWAS_CONFIG.get("mpheno", 1)
GWAS_QCOVAR = GWAS_CONFIG.get("qcovar")
GWAS_NPC = int(GWAS_CONFIG.get("npc", 10))
GWAS_HAS_COVAR = GWAS_NPC > 0 or bool(GWAS_QCOVAR)
GWAS_SPARSE_CUTOFF = GWAS_CONFIG.get("sparse_cutoff", 0.05)
GWAS_FRACTION = float(GWAS_CONFIG.get("gwas_fraction", 0.5))
GWAS_SEED = int(GWAS_CONFIG.get("seed", 42))
GWAS_PR_CONFIG = GWAS_CONFIG.get("prs", {})
GWAS_PR_ENABLED = bool(GWAS_PR_CONFIG.get("enabled", False))
GWAS_PR_SOURCE_SUBSET = GWAS_PR_CONFIG.get(
    "source_subset", GWAS_CONFIG.get("subset", "full")
)
GWAS_PR_TARGET_SUBSET = GWAS_PR_CONFIG.get("target_subset", GWAS_PR_SOURCE_SUBSET)

TRACTOR_CONFIG = GWAS_CONFIG.get("tractor", {})
TRACTOR_NUM_ANCS = int(TRACTOR_CONFIG.get("num_ancestries", 2))
TRACTOR_METHOD = TRACTOR_CONFIG.get("method", "linear")
TRACTOR_VCF_PATTERN = TRACTOR_CONFIG.get("vcf_pattern", "chr{CHR}.phased.vcf.gz")
TRACTOR_MSP_PATTERN = TRACTOR_CONFIG.get("msp_pattern", "chr{CHR}.lai.msp.tsv")
TRACTOR_DIR = TRACTOR_CONFIG.get("dir", "02-localAncestry")

GWAS_PLINK2_CONTAINER = "docker://gfanz/plink2:latest"
# Use local SIF during development/testing; push to GHCR once verified.
GWAS_TRACTOR_CONTAINER = "file://" + str((Path(workflow.basedir) / ".." / "envs" / "tractor.sif").resolve())

if GWAS_ACTIVE:
    valid_methods = {"standard", "fastgwa", "tractor"}
    unknown = set(GWAS_METHODS) - valid_methods
    if unknown:
        raise ValueError(f"gwas.methods must be in {valid_methods}, got {unknown}")
    if GWAS_FASTGWA_MODE not in {"lr", "mlm"}:
        raise ValueError("gwas.fastgwa_mode must be 'lr' or 'mlm'")
    if "tractor" in GWAS_METHODS and not uses_local_ancestry():
        raise ValueError(
            "gwas.methods includes 'tractor' but no local ancestry output is available"
            " (set ancestry.model to rfmix/gnomix)"
        )
    if GWAS_PR_ENABLED and GWAS_PR_TARGET_SUBSET != GWAS_PR_SOURCE_SUBSET:
        raise ValueError(
            "gwas.prs.target_subset must equal gwas.prs.source_subset;"
            " use prsPipeline.external for cross-ancestry transfer"
        )


def gwas_subsets():
    """Populations to test: explicit config list, else full + inferred ancestries."""
    explicit = GWAS_CONFIG.get("subsets")
    if explicit:
        return list(explicit)
    inferred = get_ancestries(None)
    return ["full"] + [a for a in inferred if a != "full"]


def _gwas_outdir(subset):
    return OUT_DIR / subset / GWAS_OUTDIR


def _gwas_base(subset):
    return OUT_DIR / subset / GWAS_BED_PREFIX


def _gwas_eigenvec(wildcards):
    template = GWAS_CONFIG.get("eigenvec")
    if template:
        return Path(str(template).format(subset=wildcards.subset))
    return OUT_DIR / wildcards.subset / "internal_pca_plink2.eigenvec"


def _tractor_input(pattern, wildcards):
    path = Path(str(pattern).format(CHR=wildcards.CHR))
    return path if path.is_absolute() else OUT_DIR / TRACTOR_DIR / path


def _vcf_prefix_template(pattern):
    """Filename prefix (with {CHR}) that extract_tracts.py prepends to outputs."""
    name = Path(str(pattern)).name
    for suffix in (".gz", ".bgz"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    if name.endswith(".vcf"):
        name = name[: -len(".vcf")]
    return name


_TRACTOR_PREFIX_TEMPLATE = _vcf_prefix_template(TRACTOR_VCF_PATTERN)


def _gwas_method_sumstats(subset, method, ancestry=None):
    base = _gwas_outdir(subset)
    if method == "standard":
        return base / "standard.sumstats.tsv"
    if method == "fastgwa":
        return base / "fastgwa.sumstats.tsv"
    if method == "tractor":
        return base / f"tractor.anc{ancestry}.sumstats.tsv"
    raise ValueError(method)


def _gwas_outputs_for_subset(subset):
    outs = []
    for method in GWAS_METHODS:
        if method == "tractor":
            outs.extend(
                _gwas_outdir(subset) / f"tractor.anc{i}.sumstats.tsv"
                for i in range(TRACTOR_NUM_ANCS)
            )
        else:
            outs.append(_gwas_method_sumstats(subset, method))
    return outs


if GWAS_ACTIVE:

    rule convertToBedForGwas:
        log:
            OUT_DIR / "logs" / "convertToBedForGwas_{subset}.log",
        container:
            GWAS_PLINK2_CONTAINER,
        conda:
            "../../envs/ancNreport.yml",
        envmodules:
            *([config.get("plink_module")] if config.get("plink_module") else []),
        threads: 8,
        resources:
            nodes=1,
            mem_mb=32000,
            runtime=60,
        input:
            pgen=lambda w: Path(str(_gwas_base(w.subset)) + ".pgen"),
            pvar=lambda w: Path(str(_gwas_base(w.subset)) + ".pvar"),
            psam=lambda w: Path(str(_gwas_base(w.subset)) + ".psam"),
        output:
            bed=OUT_DIR / "{subset}" / (GWAS_BED_PREFIX + ".bed"),
            bim=OUT_DIR / "{subset}" / (GWAS_BED_PREFIX + ".bim"),
            fam=OUT_DIR / "{subset}" / (GWAS_BED_PREFIX + ".fam"),
        params:
            in_prefix=lambda w, input: str(input.pgen)[:-5],
            out_prefix=lambda w, output: str(output.bed)[:-4],
        shell:
            """
            plink2 --pfile {params.in_prefix} --make-bed --out {params.out_prefix} --threads {threads}
            """

    rule gwasSampleSplit:
        log:
            OUT_DIR / "logs" / "gwasSampleSplit_{subset}.log",
        conda:
            "../../envs/ancNreport.yml",
        threads: 1,
        resources:
            nodes=1,
            mem_mb=4000,
            runtime=10,
        input:
            fam=OUT_DIR / "{subset}" / (GWAS_BED_PREFIX + ".fam"),
        output:
            gwas_keep=_gwas_outdir("{subset}") / "samples.gwas.keep",
            prs_keep=_gwas_outdir("{subset}") / "samples.prs.keep",
        params:
            fraction=GWAS_FRACTION,
            seed=GWAS_SEED,
            split_flag="--split" if GWAS_PR_ENABLED else "",
        shell:
            """
            python {SCRIPTS_DIR}/select_gwas_samples.py \
                --fam {input.fam} \
                --gwas-keep {output.gwas_keep} \
                --prs-keep {output.prs_keep} \
                --fraction {params.fraction} \
                --seed {params.seed} \
                {params.split_flag}
            """

    rule prepareGwasInputs:
        log:
            OUT_DIR / "logs" / "prepareGwasInputs_{subset}.log",
        conda:
            "../../envs/ancNreport.yml",
        threads: 1,
        resources:
            nodes=1,
            mem_mb=4000,
            runtime=10,
        input:
            pheno=GWAS_PHENO,
            eigenvec=lambda w: (
                [_gwas_eigenvec(w)] if GWAS_NPC > 0 else []
            ),
            qcovar=lambda w: [GWAS_QCOVAR] if GWAS_QCOVAR else [],
            gwas_keep=_gwas_outdir("{subset}") / "samples.gwas.keep",
            prs_keep=_gwas_outdir("{subset}") / "samples.prs.keep",
        output:
            gwas_pheno=_gwas_outdir("{subset}") / "gwas_pheno.tsv",
            prs_pheno=_gwas_outdir("{subset}") / "prs_pheno.tsv",
            gwas_covar=_gwas_outdir("{subset}") / "gwas_covar.tsv",
            prs_covar=_gwas_outdir("{subset}") / "prs_covar.tsv",
            prs_pcs=_gwas_outdir("{subset}") / "prs_pcs.tsv",
            tractor_pheno=_gwas_outdir("{subset}") / "tractor_pheno.tsv",
            covar_columns=_gwas_outdir("{subset}") / "gwas_covar.tsv.columns",
        params:
            mpheno=GWAS_MPHENO,
            npc=GWAS_NPC,
            eigenvec_flag=lambda w, input: (
                f"--eigenvec {input.eigenvec[0]}" if input.eigenvec else ""
            ),
            qcovar_flag=(
                f"--qcovar {GWAS_QCOVAR}" if GWAS_QCOVAR else ""
            ),
        shell:
            """
            python {SCRIPTS_DIR}/prepare_gwas_inputs.py \
                --pheno {input.pheno} \
                --mpheno {params.mpheno} \
                --npc {params.npc} \
                {params.eigenvec_flag} \
                {params.qcovar_flag} \
                --gwas-keep {input.gwas_keep} \
                --prs-keep {input.prs_keep} \
                --gwas-pheno {output.gwas_pheno} \
                --prs-pheno {output.prs_pheno} \
                --gwas-covar {output.gwas_covar} \
                --prs-covar {output.prs_covar} \
                --prs-pcs {output.prs_pcs} \
                --tractor-pheno {output.tractor_pheno}
            """

    if "standard" in GWAS_METHODS:

        rule runStandardGwas:
            log:
                OUT_DIR / "logs" / "runStandardGwas_{subset}.log",
            container:
                GWAS_PLINK2_CONTAINER,
            conda:
                "../../envs/ancNreport.yml",
            envmodules:
                *([config.get("plink_module")] if config.get("plink_module") else []),
            threads: 8,
            resources:
                nodes=1,
                mem_mb=32000,
                runtime=480,
            input:
                bed=OUT_DIR / "{subset}" / (GWAS_BED_PREFIX + ".bed"),
                bim=OUT_DIR / "{subset}" / (GWAS_BED_PREFIX + ".bim"),
                fam=OUT_DIR / "{subset}" / (GWAS_BED_PREFIX + ".fam"),
                keep=_gwas_outdir("{subset}") / "samples.gwas.keep",
                pheno=_gwas_outdir("{subset}") / "gwas_pheno.tsv",
                covar=_gwas_outdir("{subset}") / "gwas_covar.tsv",
            output:
                raw=_gwas_outdir("{subset}") / "standard.PHENO.glm.linear",
                sumstats=_gwas_outdir("{subset}") / "standard.sumstats.tsv",
            params:
                bed_prefix=lambda w: str(_gwas_base(w.subset)),
                out_prefix=lambda w: str(_gwas_outdir(w.subset) / "standard"),
                covar_flag=(
                    lambda w, input: f"--covar {input.covar}"
                    if GWAS_HAS_COVAR
                    else ""
                ),
            shell:
                """
                plink2 --bfile {params.bed_prefix} \
                    --keep {input.keep} \
                    --pheno {input.pheno} \
                    --pheno-name PHENO \
                    {params.covar_flag} \
                    --glm hide-covar allow-no-covars \
                    --threads {threads} \
                    --out {params.out_prefix}
                python {SCRIPTS_DIR}/format_gwas_sumstats.py \
                    --method plink2 \
                    --input {output.raw} \
                    --output {output.sumstats}
                """

    if "fastgwa" in GWAS_METHODS:

        if GWAS_FASTGWA_MODE == "mlm":

            rule makeGwasSparseGrm:
                log:
                    OUT_DIR / "logs" / "makeGwasSparseGrm_{subset}.log",
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
                    grm_bin=lambda w: (
                        OUT_DIR / w.subset / (GWAS_GRM_PREFIX + ".grm.bin")
                    ),
                    keep=_gwas_outdir("{subset}") / "samples.gwas.keep",
                output:
                    sparse=_gwas_outdir("{subset}") / "gwas_sparsegrm.grm.sp",
                params:
                    grm_prefix=lambda w: str(OUT_DIR / w.subset / GWAS_GRM_PREFIX),
                    out_prefix=lambda w, output: str(output.sparse)[: -len(".grm.sp")],
                    cutoff=GWAS_SPARSE_CUTOFF,
                shell:
                    """
                    gcta --grm {params.grm_prefix} \
                        --keep {input.keep} \
                        --make-bK-sparse {params.cutoff} \
                        --out {params.out_prefix} \
                        --thread-num {threads}
                    """

        rule runFastGwa:
            log:
                OUT_DIR / "logs" / "runFastGwa_{subset}.log",
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
                bed=OUT_DIR / "{subset}" / (GWAS_BED_PREFIX + ".bed"),
                bim=OUT_DIR / "{subset}" / (GWAS_BED_PREFIX + ".bim"),
                fam=OUT_DIR / "{subset}" / (GWAS_BED_PREFIX + ".fam"),
                keep=_gwas_outdir("{subset}") / "samples.gwas.keep",
                pheno=_gwas_outdir("{subset}") / "gwas_pheno.tsv",
                covar=_gwas_outdir("{subset}") / "gwas_covar.tsv",
                sparse=(
                    _gwas_outdir("{subset}") / "gwas_sparsegrm.grm.sp"
                    if GWAS_FASTGWA_MODE == "mlm"
                    else []
                ),
            output:
                raw=_gwas_outdir("{subset}") / "fastgwa.fastGWA",
                sumstats=_gwas_outdir("{subset}") / "fastgwa.sumstats.tsv",
            params:
                bed_prefix=lambda w: str(_gwas_base(w.subset)),
                out_prefix=lambda w: str(_gwas_outdir(w.subset) / "fastgwa"),
                grm_flag=(
                    lambda w, input: (
                        f"--grm-sparse {str(input.sparse)[: -len('.grm.sp')]}"
                        " --fastGWA-mlm"
                    )
                    if GWAS_FASTGWA_MODE == "mlm"
                    else "--fastGWA-lr"
                ),
                qcovar_flag=(
                    lambda w, input: f"--qcovar {input.covar}"
                    if GWAS_HAS_COVAR
                    else ""
                ),
            shell:
                """
                gcta --bfile {params.bed_prefix} \
                    --keep {input.keep} \
                    --pheno {input.pheno} \
                    --mpheno 1 \
                    {params.qcovar_flag} \
                    {params.grm_flag} \
                    --out {params.out_prefix}
                python {SCRIPTS_DIR}/format_gwas_sumstats.py \
                    --method gcta \
                    --input {output.raw} \
                    --output {output.sumstats}
                """

    if "tractor" in GWAS_METHODS:

        rule extractGwasTracts:
            log:
                OUT_DIR / "logs" / "extractGwasTracts_{subset}_{CHR}.log",
            container:
                GWAS_TRACTOR_CONTAINER,
            conda:
                "../../envs/tractor.yml",
            threads: 4,
            resources:
                nodes=1,
                mem_mb=16000,
                runtime=720,
            wildcard_constraints:
                CHR=CHR_REGEX,
            input:
                vcf=lambda w: _tractor_input(TRACTOR_VCF_PATTERN, w),
                msp=lambda w: _tractor_input(TRACTOR_MSP_PATTERN, w),
            output:
                dosage=expand(
                    _gwas_outdir("{subset}")
                    / "tractor"
                    / f"{_TRACTOR_PREFIX_TEMPLATE}_anc{{anc}}.dosage.txt",
                    anc=range(TRACTOR_NUM_ANCS),
                    subset=["{subset}"],
                    CHR=["{CHR}"],
                ),
                hapcount=expand(
                    _gwas_outdir("{subset}")
                    / "tractor"
                    / f"{_TRACTOR_PREFIX_TEMPLATE}_anc{{anc}}.hapcount.txt",
                    anc=range(TRACTOR_NUM_ANCS),
                    subset=["{subset}"],
                    CHR=["{CHR}"],
                ),
            params:
                num_ancs=TRACTOR_NUM_ANCS,
                out_dir=lambda w: _gwas_outdir(w.subset) / "tractor",
                subset="{subset}",
            shell:
                """
                python /opt/Tractor/scripts/extract_tracts.py \
                    --vcf {input.vcf} \
                    --msp {input.msp} \
                    --num-ancs {params.num_ancs} \
                    --output-dir {params.out_dir}
                """

        rule mergeGwasTracts:
            log:
                OUT_DIR / "logs" / "mergeGwasTracts_{subset}.log",
            conda:
                "../../envs/tractor.yml",
            threads: 1,
            resources:
                nodes=1,
                mem_mb=8000,
                runtime=60,
            input:
                dosage=expand(
                    _gwas_outdir("{subset}")
                    / "tractor"
                    / f"{_TRACTOR_PREFIX_TEMPLATE}_anc{{anc}}.dosage.txt",
                    anc=range(TRACTOR_NUM_ANCS),
                    subset=["{subset}"],
                    CHR=LOCAL_ANCESTRY_CHROMOSOMES,
                ),
                hapcount=expand(
                    _gwas_outdir("{subset}")
                    / "tractor"
                    / f"{_TRACTOR_PREFIX_TEMPLATE}_anc{{anc}}.hapcount.txt",
                    anc=range(TRACTOR_NUM_ANCS),
                    subset=["{subset}"],
                    CHR=LOCAL_ANCESTRY_CHROMOSOMES,
                ),
            output:
                dosage=expand(
                    _gwas_outdir("{subset}") / "tractor" / "merged" / "all.anc{anc}.dosage.txt",
                    anc=range(TRACTOR_NUM_ANCS),
                    subset=["{subset}"],
                ),
                hapcount=expand(
                    _gwas_outdir("{subset}") / "tractor" / "merged" / "all.anc{anc}.hapcount.txt",
                    anc=range(TRACTOR_NUM_ANCS),
                    subset=["{subset}"],
                ),
            params:
                num_ancs=TRACTOR_NUM_ANCS,
                src_dir=lambda w: _gwas_outdir(w.subset) / "tractor",
                out_dir=lambda w: _gwas_outdir(w.subset) / "tractor" / "merged",
            shell:
                """
                mkdir -p {params.out_dir}
                for anc in $(seq 0 $(({params.num_ancs} - 1))); do
                    awk 'FNR == 1 && NR != 1 {{next}} {{print}}' \
                        {params.src_dir}/*.anc${{anc}}.dosage.txt \
                        > {params.out_dir}/all.anc${{anc}}.dosage.txt
                    awk 'FNR == 1 && NR != 1 {{next}} {{print}}' \
                        {params.src_dir}/*.anc${{anc}}.hapcount.txt \
                        > {params.out_dir}/all.anc${{anc}}.hapcount.txt
                done
                """

        rule runTractor:
            log:
                OUT_DIR / "logs" / "runTractor_{subset}.log",
            container:
                GWAS_TRACTOR_CONTAINER,
            conda:
                "../../envs/tractor.yml",
            threads: 8,
            resources:
                nodes=1,
                mem_mb=32000,
                runtime=1440,
            input:
                dosage=expand(
                    _gwas_outdir("{subset}") / "tractor" / "merged" / "all.anc{anc}.dosage.txt",
                    anc=range(TRACTOR_NUM_ANCS),
                    subset=["{subset}"],
                ),
                hapcount=expand(
                    _gwas_outdir("{subset}") / "tractor" / "merged" / "all.anc{anc}.hapcount.txt",
                    anc=range(TRACTOR_NUM_ANCS),
                    subset=["{subset}"],
                ),
                pheno=_gwas_outdir("{subset}") / "tractor_pheno.tsv",
                covar_columns=_gwas_outdir("{subset}") / "gwas_covar.tsv.columns",
            output:
                raw=_gwas_outdir("{subset}") / "tractor.raw.tsv",
                sumstats=expand(
                    _gwas_outdir("{subset}") / "tractor.anc{anc}.sumstats.tsv",
                    anc=range(TRACTOR_NUM_ANCS),
                    subset=["{subset}"],
                ),
            params:
                hapdose=lambda w: str(_gwas_outdir(w.subset) / "tractor" / "merged" / "all"),
                out_prefix=lambda w: str(_gwas_outdir(w.subset) / "tractor.raw"),
                sumstats_prefix=lambda w: str(
                    _gwas_outdir(w.subset) / "tractor.anc"
                ),
                num_ancs=TRACTOR_NUM_ANCS,
                method=TRACTOR_METHOD,
                covar_flag=(
                    lambda w, input: (
                        '--covarcollist "$(paste -sd, '
                        + str(input.covar_columns)
                        + ')"'
                    )
                    if GWAS_HAS_COVAR
                    else ""
                ),
            shell:
                """
                Rscript /opt/Tractor/scripts/run_tractor.R \
                    --hapdose {params.hapdose} \
                    --phenofile {input.pheno} \
                    --sampleidcol IID \
                    --phenocol PHENO \
                    --method {params.method} \
                    --nthreads {threads} \
                    {params.covar_flag} \
                    --output {output.raw}
                for anc in $(seq 0 $(({params.num_ancs} - 1))); do
                    python {SCRIPTS_DIR}/format_gwas_sumstats.py \
                        --method tractor \
                        --ancestry ${{anc}} \
                        --input {output.raw} \
                        --output {params.sumstats_prefix}${{anc}}.sumstats.tsv
                done
                """

    if GWAS_PR_ENABLED:
        if "standard" in GWAS_METHODS:
            _gwas_pr_default_method = "standard"
        elif "fastgwa" in GWAS_METHODS:
            _gwas_pr_default_method = "fastgwa"
        else:
            _gwas_pr_default_method = "tractor"
        GWAS_PR_METHOD = str(
            GWAS_PR_CONFIG.get("method", _gwas_pr_default_method)
        ).lower()
        if GWAS_PR_METHOD not in GWAS_METHODS:
            raise ValueError(
                f"gwas.prs.method '{GWAS_PR_METHOD}' is not in gwas.methods {GWAS_METHODS}"
            )
        GWAS_PR_TRACTOR_ANC = int(GWAS_PR_CONFIG.get("tractor_ancestry", 0))
        if GWAS_PR_METHOD == "tractor" and GWAS_PR_TRACTOR_ANC >= TRACTOR_NUM_ANCS:
            raise ValueError(
                f"gwas.prs.tractor_ancestry must be < tractor.num_ancestries"
                f" ({TRACTOR_NUM_ANCS})"
            )
        _gwas_pr_base = _gwas_outdir(GWAS_PR_SOURCE_SUBSET)

        rule stageGwasForPrs:
            log:
                OUT_DIR / "logs" / "stageGwasForPrs.log",
            container:
                GWAS_PLINK2_CONTAINER,
            conda:
                "../../envs/ancNreport.yml",
            envmodules:
                *([config.get("plink_module")] if config.get("plink_module") else []),
            threads: 4,
            resources:
                nodes=1,
                mem_mb=16000,
                runtime=120,
            input:
                bed=OUT_DIR / GWAS_PR_TARGET_SUBSET / (GWAS_BED_PREFIX + ".bed"),
                bim=OUT_DIR / GWAS_PR_TARGET_SUBSET / (GWAS_BED_PREFIX + ".bim"),
                fam=OUT_DIR / GWAS_PR_TARGET_SUBSET / (GWAS_BED_PREFIX + ".fam"),
                keep=_gwas_outdir(GWAS_PR_TARGET_SUBSET) / "samples.prs.keep",
                pheno=_gwas_outdir(GWAS_PR_TARGET_SUBSET) / "prs_pheno.tsv",
                pcs=_gwas_outdir(GWAS_PR_TARGET_SUBSET) / "prs_pcs.tsv",
                sumstats=_gwas_method_sumstats(
                    GWAS_PR_SOURCE_SUBSET, GWAS_PR_METHOD, GWAS_PR_TRACTOR_ANC
                ),
            output:
                bed=PRS_OUT_DIR / "gwas_staged" / "target.bed",
                bim=PRS_OUT_DIR / "gwas_staged" / "target.bim",
                fam=PRS_OUT_DIR / "gwas_staged" / "target.fam",
                pheno=PRS_OUT_DIR / "gwas_staged" / "target_pheno.tsv",
                pcs=PRS_OUT_DIR / "gwas_staged" / "target_pcs.tsv",
                sumstats=PRS_OUT_DIR / "gwas_staged" / "target_sumstats.tsv",
                ready=PRS_OUT_DIR / "gwas_staged" / "staged.ready",
            params:
                bed_prefix=lambda w: str(
                    OUT_DIR / GWAS_PR_TARGET_SUBSET / GWAS_BED_PREFIX
                ),
                out_prefix=lambda w: str(PRS_OUT_DIR / "gwas_staged" / "target"),
            shell:
                """
                plink2 --bfile {params.bed_prefix} \
                    --keep {input.keep} \
                    --make-bed \
                    --threads {threads} \
                    --out {params.out_prefix}
                cp {input.pheno} {output.pheno}
                cp {input.pcs} {output.pcs}
                cp {input.sumstats} {output.sumstats}
                touch {output.ready}
                """


if GWAS_ACTIVE:

    rule run_gwas:
        input:
            lambda wildcards: sum(
                (_gwas_outputs_for_subset(s) for s in gwas_subsets()), []
            ),
