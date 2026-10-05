import shlex
import shutil

PRS_METHODS_CONFIG = config.get("prsMethods", {})

# Absorbed prs_pipeline engine = vendored copy under workflow/scripts/prs_pipeline
# (see VENDORED.md). Resolved/presence-checked in preparePRSInputs.smk.

# SIF container paths used by upstream prs_pipeline method rules. Resolve from
# config under prsMethods.containers (keys: prsv2, singleprshelper); default to a
# Phase-7 landing spot under OUT_DIR/containers (the SIFs are not vendored). The
# primary image is the GDCGenomicsQC-owned `prs` build (envs/prs.def), which bakes
# the vendored engine + pinned method repos into /opt/gdcgenomicsqc. Pull:
#   apptainer pull oras://ghcr.io/mainsqu33ze/gdcgenomicsqc/prs:latest
#   apptainer pull oras://ghcr.io/mainsqu33ze/gdcgenomicsqc/singleprshelper:latest
PRS_CONTAINERS = {
    "prsv2": PRS_METHODS_CONFIG.get("containers", {}).get(
        "prsv2", str(OUT_DIR / "containers" / "prs_latest.sif")
    ),
    "singleprshelper": PRS_METHODS_CONFIG.get("containers", {}).get(
        "singleprshelper", str(OUT_DIR / "containers" / "singleprshelper_latest.sif")
    ),
}
PRS_RESOURCE_DIR = Path(
    PRS_METHODS_CONFIG.get(
        "resource_dir",
        str(OUT_DIR.parent / "prs_resources"),
    )
)
PRS_METHOD_RUN_DIR = PRS_OUT_DIR / "method_runs"
PRS_DOWNLOAD_SOFTWARE_FLAG = (
    "--download-software" if PRS_METHODS_CONFIG.get("download_software", False) else ""
)

# --- Phase 2: native single-ancestry method plumbing (vendored engine) ---
PRS_SRC = PRS_PIPELINE_PATH / "src"
PRS_STUDY_PLINK = PRS_OUT_DIR / "anc1_plink_files" / f"{PRS_ANC1}_simulation_study_sample"
PRS_STUDY_RDS = PRS_STUDY_PLINK.with_suffix(".rds")
PRS_STUDY_BK = PRS_STUDY_PLINK.with_suffix(".bk")
# Shared ALIGNED summary stats (monolithic step 1 output) + study-sample phenotype
# (monolithic step 2 output). Both feed CT / LDpred2 / lassosum2 / PRSice2.
PRS_ALIGNED_SUMSTATS = PRS_OUT_DIR / "gwas" / "CT_PRSice2_summary_stat_file.txt"
PRS_STUDY_PHENO_FILE = PRS_OUT_DIR / "gwas" / "study_sample_pheno.txt"
# Shared per-ancestry LD-matrix dir; per-method config (single_ldpred2.ld_matrix_dir
# etc.) overrides the shared default with an externally-supplied matrix.
PRS_LD_MATRIX_DIR = Path(
    PRS_CONFIG.get("ld_matrix_dir") or str(PRS_OUT_DIR / "ld_matrix" / PRS_ANC1)
)
PRS_N_TOTAL_GWAS = PRS_CONFIG.get("n_total_gwas", 31968)
PRS_AFREQ_FILE = PRS_CONFIG.get("afreq_file", "")
PRS_PLINK_BIN = PRS_CONFIG.get("path_plink", "plink")
PRS_PLINK2_BIN = PRS_CONFIG.get("path_plink2", "") or PRS_METHODS_CONFIG.get("plink2", "") or "plink2"
# Single-ancestry PRS-CS (PRScsx.py single-pop). Defaults to PRS-CSx code + LD
# reference provisioned under the resource dir (download_prs_resources.sh);
# override via prsMethods.single_prscs.path_code / ld_ref_dir / seed.
PRS_PRSCS_PATH_CODE = (
    PRS_METHODS_CONFIG.get("single_prscs", {}).get("path_code")
    or str(PRS_RESOURCE_DIR / "software" / "PRScsx")
)
PRS_PRSCS_REF_DIR = (
    PRS_METHODS_CONFIG.get("single_prscs", {}).get("ld_ref_dir")
    or str(PRS_RESOURCE_DIR / "ld" / "prs_csx" / "ref")
)
PRS_PRSCS_SEED = PRS_METHODS_CONFIG.get("single_prscs", {}).get("seed", 42)
PRS_PRSCS_PATH_PYTHON = PRS_METHODS_CONFIG.get("single_prscs", {}).get("path_python", None)
# Joint multi-ancestry PRS-CSx (PRScsx.py in multi-population mode). Reuses the
# single_prscs code/LD-reference/seed defaults unless overridden via
# prsMethods.multi_prscsx.path_code / ld_ref_dir / seed / path_python.
PRS_PRSCSX_PATH_CODE = (
    PRS_METHODS_CONFIG.get("multi_prscsx", {}).get("path_code")
    or PRS_PRSCS_PATH_CODE
)
PRS_PRSCSX_REF_DIR = (
    PRS_METHODS_CONFIG.get("multi_prscsx", {}).get("ld_ref_dir")
    or PRS_PRSCS_REF_DIR
)
PRS_PRSCSX_SEED = PRS_METHODS_CONFIG.get("multi_prscsx", {}).get("seed", PRS_PRSCS_SEED)
PRS_PRSCSX_PATH_PYTHON = PRS_METHODS_CONFIG.get(
    "multi_prscsx", {}
).get("path_python", PRS_PRSCS_PATH_PYTHON)

# --- Phase 3: test evaluation (score_test.sh) ---
PRS_TEST_SAMPLE = PRS_CONFIG.get("test_sample", "")
PRS_TEST_PCA_FILE = PRS_CONFIG.get("test_pca_eigenvec_file", "")
PRS_TEST_EVAL_DIR = PRS_METHOD_RUN_DIR / "test_evaluation"
# Staging dir where score_test.sh's upstream layout is bridged to the native
# per-method layout (see scoreTestPRS); score_test writes into
# <stage>/test_evaluation, which is then moved up to PRS_TEST_EVAL_DIR.
PRS_TEST_STAGE_DIR = PRS_METHOD_RUN_DIR / "test_evaluation" / "_stage"
PRS_TEST_INPUTS = (
    {
        "test_bed": Path(PRS_TEST_SAMPLE).with_suffix(".bed"),
        "test_bim": Path(PRS_TEST_SAMPLE).with_suffix(".bim"),
        "test_fam": Path(PRS_TEST_SAMPLE).with_suffix(".fam"),
    }
    if PRS_TEST_SAMPLE
    else {}
)


def prs_method_ld_matrix_dir(method):
    return PRS_METHODS_CONFIG.get(method, {}).get("ld_matrix_dir") or str(PRS_LD_MATRIX_DIR)


def prs_method_command(method):
    return PRS_METHODS_CONFIG.get(method, {}).get("command", "")


def prs_method_command_quoted(method):
    return shlex.quote(prs_method_command(method))


def prs_method_extra_args(method):
    method_config = PRS_METHODS_CONFIG.get(method, {})
    args = []
    for key, flag in (
        ("ld_ref_dir", "--ld-ref-dir"),
        ("ld_ref_prefix", "--ld-ref-prefix"),
        ("ld_matrix_dir", "--ld-matrix-dir"),
        ("software_dir", "--software-dir"),
    ):
        value = method_config.get(key)
        if value:
            args.extend([flag, shlex.quote(str(value))])
    return " ".join(args)


rule preparePRSMethodResources:
    log:
        OUT_DIR / "logs" / "preparePRSMethodResources.log",
    container:
        "oras://ghcr.io/mainsqu33ze/gdcgenomicsqc/prs:latest"
    conda:
        "../../envs/prs.yml"
    threads: 1
    resources:
        nodes=1,
        mem_mb=4000,
        runtime=60,
    input:
        prs_pipeline_ready=rules.checkPRSPipelinePath.output.checked,
    output:
        ready=PRS_RESOURCE_DIR / "resources.ready",
    params:
        resource_dir=PRS_RESOURCE_DIR,
        script=Path(workflow.basedir) / "scripts" / "download_prs_resources.sh",
        prs_pipeline_dir=PRS_PIPELINE_PATH,
        prs_pipeline_ref=PRS_PIPELINE_PIN,
        prs_pipeline_sif=PRS_CONTAINERS["prsv2"],
        prs_helper_sif=PRS_CONTAINERS["singleprshelper"],
        prscsx_ref=PRS_CONFIG.get(
            "path_ref_dir",
            PRS_METHODS_CONFIG.get("multi_prscsx", {}).get("ld_ref_dir", ""),
        ),
        plink2=PRS_CONFIG.get("path_plink2", PRS_METHODS_CONFIG.get("plink2", "")),
        download_software=PRS_DOWNLOAD_SOFTWARE_FLAG,
    shell:
        """
        bash {params.script} \
            --resource-dir {params.resource_dir} \
            --prscsx-ref-dir {params.prscsx_ref} \
            --plink2 {params.plink2} \
            --prs-pipeline-dir {params.prs_pipeline_dir} \
            --prs-pipeline-ref {params.prs_pipeline_ref} \
            --prs-pipeline-sif {params.prs_pipeline_sif} \
            --prs-helper-sif {params.prs_helper_sif} \
            {params.download_software} \
            > {log} 2>&1
        """


rule alignSumstatsForPRS:
    log:
        OUT_DIR / "logs" / "alignSumstatsForPRS.log",
    container:
        "oras://ghcr.io/mainsqu33ze/gdcgenomicsqc/prs:latest"
    conda:
        "../../envs/prs.yml"
    threads: 1
    resources:
        nodes=1,
        mem_mb=8000,
        runtime=60,
    input:
        ss=rules.preparePRSInputs.output.target_single_sumstats,
        bim=rules.preparePRSInputs.output.study_bim,
    output:
        aligned=PRS_ALIGNED_SUMSTATS,
    params:
        n_total=PRS_N_TOTAL_GWAS,
        script=PRS_SRC / "prepare_sumstats.R",
    shell:
        """
        Rscript {params.script} \\
            --input {input.ss} \\
            --bim {input.bim} \\
            --n_total {params.n_total} \\
            --output {output.aligned} > {log} 2>&1
        """


rule makeStudyPhenoFile:
    log:
        OUT_DIR / "logs" / "makeStudyPhenoFile.log",
    threads: 1
    resources:
        nodes=1,
        mem_mb=2000,
        runtime=15,
    input:
        fam=rules.preparePRSInputs.output.study_fam,
    output:
        pheno=PRS_STUDY_PHENO_FILE,
    shell:
        """
        awk 'BEGIN {{print "FID\\tIID\\tphenotype"}} {{print $1, $2, $6}}' OFS="\\t" {input.fam} > {output.pheno}
        """


rule convertStudyBedToRDS:
    log:
        OUT_DIR / "logs" / "convertStudyBedToRDS.log",
    container:
        "oras://ghcr.io/mainsqu33ze/gdcgenomicsqc/prs:latest"
    conda:
        "../../envs/prs.yml"
    threads: 1
    resources:
        nodes=1,
        mem_mb=8000,
        runtime=120,
    input:
        bed=rules.preparePRSInputs.output.study_bed,
    output:
        rds=PRS_STUDY_RDS,
        bk=PRS_STUDY_BK,
    shell:
        """
        set -euo pipefail
        if [[ ! -f "{output.rds}" ]]; then
            rm -f "{output.rds}" "{output.bk}"
        fi
        Rscript -e 'library(bigsnpr); snp_readBed("{input.bed}")' > {log} 2>&1
        """


rule generateLDMatrix:
    log:
        OUT_DIR / "logs" / "generateLDMatrix.log",
    container:
        "oras://ghcr.io/mainsqu33ze/gdcgenomicsqc/prs:latest"
    conda:
        "../../envs/prs.yml"
    threads: 4
    resources:
        nodes=1,
        mem_mb=32000,
        runtime=600,
    input:
        bed=rules.preparePRSInputs.output.study_bed,
        rds=rules.convertStudyBedToRDS.output.rds,
    output:
        map=PRS_LD_MATRIX_DIR / "map.rds",
        g_idx=PRS_LD_MATRIX_DIR / "g_idx.rds",
    params:
        script=PRS_SRC / "generate_ld_matrix.R",
        out_dir=PRS_LD_MATRIX_DIR,
    shell:
        """
        mkdir -p {params.out_dir}
        Rscript {params.script} \\
            --anc_bed {input.bed} \\
            --out {params.out_dir} \\
            --ncores {threads} > {log} 2>&1
        """


rule runSingleAncestryCT:
    log:
        OUT_DIR / "logs" / "runSingleAncestryCT.log",
    container:
        "oras://ghcr.io/mainsqu33ze/gdcgenomicsqc/prs:latest"
    conda:
        "../../envs/prs.yml"
    threads: 4
    resources:
        nodes=1,
        mem_mb=16000,
        runtime=240,
    input:
        resources=rules.preparePRSMethodResources.output.ready,
        ss=rules.alignSumstatsForPRS.output.aligned,
        pheno=rules.makeStudyPhenoFile.output.pheno,
        study_bed=rules.preparePRSInputs.output.study_bed,
        study_bim=rules.preparePRSInputs.output.study_bim,
        study_fam=rules.preparePRSInputs.output.study_fam,
    output:
        results=PRS_METHOD_RUN_DIR / "single_ct" / "CT" / "CT_prs_results.txt",
        done=PRS_METHOD_RUN_DIR / "single_ct.done",
    params:
        study_prefix=PRS_STUDY_PLINK,
        out_base=PRS_METHOD_RUN_DIR / "single_ct",
        plink=PRS_PLINK_BIN,
        pca=PRS_CONFIG.get("gwas_pca_eigenvec_file", ""),
        script=PRS_SRC / "run_CT.sh",
    shell:
        """
        set -euo pipefail

        mkdir -p {params.out_base}/CT/temp

        CONFIG={params.out_base}/CT/temp/CT_temp_config.txt
        echo "study_sample={params.study_prefix}" > "$CONFIG"
        echo "sum_stats_file={input.ss}" >> "$CONFIG"
        echo "phenotype_info_file={input.pheno}" >> "$CONFIG"
        echo "output_path={params.out_base}" >> "$CONFIG"
        echo "path_prs_pipeline={PRS_PIPELINE_PATH}" >> "$CONFIG"
        if [[ -n "{params.pca}" ]]; then
            echo "gwas_pca_eigenvec_file={params.pca}" >> "$CONFIG"
        fi

        PLINK_BIN="{params.plink}"
        PLINK_DIR="$(dirname "$PLINK_BIN")"
        if [[ "$PLINK_BIN" != "plink" && "$PLINK_DIR" != "." ]]; then
            export PATH="$PLINK_DIR:$PATH"
        fi

        bash {params.script} --c "$CONFIG" > {log} 2>&1

        touch {output.done}
        """


rule runSingleAncestryLDpred2:
    log:
        OUT_DIR / "logs" / "runSingleAncestryLDpred2.log",
    container:
        "oras://ghcr.io/mainsqu33ze/gdcgenomicsqc/prs:latest"
    conda:
        "../../envs/prs.yml"
    threads: 4
    resources:
        nodes=1,
        mem_mb=64000,
        runtime=720,
    input:
        resources=rules.preparePRSMethodResources.output.ready,
        ss=rules.alignSumstatsForPRS.output.aligned,
        pheno=rules.makeStudyPhenoFile.output.pheno,
        rds=rules.convertStudyBedToRDS.output.rds,
        bk=rules.convertStudyBedToRDS.output.bk,
        study_bed=rules.preparePRSInputs.output.study_bed,
        study_bim=rules.preparePRSInputs.output.study_bim,
        study_fam=rules.preparePRSInputs.output.study_fam,
        ld_map=lambda wildcards: [
            rules.generateLDMatrix.output.map
            if prs_method_ld_matrix_dir("single_ldpred2") == str(PRS_LD_MATRIX_DIR)
            else []
        ][0],
    output:
        scores=PRS_METHOD_RUN_DIR / "single_ldpred2" / "prs_method_individual_scores.txt",
        performance=PRS_METHOD_RUN_DIR / "single_ldpred2" / "prs_method_performance.csv",
        inf_weights=PRS_METHOD_RUN_DIR / "single_ldpred2" / "prs_method_inf_weights.txt",
        grid_weights=PRS_METHOD_RUN_DIR / "single_ldpred2" / "prs_method_grid_weights.txt",
        plot=PRS_METHOD_RUN_DIR / "single_ldpred2" / "prs_method_grid_plot.png",
        done=PRS_METHOD_RUN_DIR / "single_ldpred2.done",
    params:
        out_base=PRS_METHOD_RUN_DIR / "single_ldpred2",
        out_prefix=PRS_METHOD_RUN_DIR / "single_ldpred2" / "prs_method",
        ld_matrix_dir=lambda wildcards: prs_method_ld_matrix_dir("single_ldpred2"),
        afreq=PRS_AFREQ_FILE,
        script=PRS_SRC / "run_LDpred2.R",
    shell:
        """
        set -euo pipefail

        mkdir -p {params.out_base}
        cd {params.out_base}

        LDpred2_args="--rds {input.rds} --ss {input.ss} --bim {input.study_bim} --out {params.out_prefix} --ncores {threads} --pheno {input.pheno} --ld-matrix-dir {params.ld_matrix_dir}"
        if [[ -n "{params.afreq}" ]]; then
            LDpred2_args="$LDpred2_args --afreq {params.afreq}"
        fi

        Rscript {params.script} $LDpred2_args > {log} 2>&1

        touch {output.done}
        """


rule runSingleAncestryLassosum2:
    log:
        OUT_DIR / "logs" / "runSingleAncestryLassosum2.log",
    container:
        "oras://ghcr.io/mainsqu33ze/gdcgenomicsqc/prs:latest"
    conda:
        "../../envs/prs.yml"
    threads: 4
    resources:
        nodes=1,
        mem_mb=64000,
        runtime=720,
    input:
        resources=rules.preparePRSMethodResources.output.ready,
        ss=rules.alignSumstatsForPRS.output.aligned,
        pheno=rules.makeStudyPhenoFile.output.pheno,
        rds=rules.convertStudyBedToRDS.output.rds,
        bk=rules.convertStudyBedToRDS.output.bk,
        study_bed=rules.preparePRSInputs.output.study_bed,
        study_bim=rules.preparePRSInputs.output.study_bim,
        study_fam=rules.preparePRSInputs.output.study_fam,
        ld_map=lambda wildcards: [
            rules.generateLDMatrix.output.map
            if prs_method_ld_matrix_dir("single_lassosum2") == str(PRS_LD_MATRIX_DIR)
            else []
        ][0],
    output:
        predictions=PRS_METHOD_RUN_DIR / "single_lassosum2" / "prs_method_full_predictions.csv",
        best_prs=PRS_METHOD_RUN_DIR / "single_lassosum2" / "prs_method_final_best_prs.csv",
        weights=PRS_METHOD_RUN_DIR / "single_lassosum2" / "prs_method_weights.txt",
        grid_params=PRS_METHOD_RUN_DIR / "single_lassosum2" / "prs_method_grid_params.csv",
        final_res=PRS_METHOD_RUN_DIR / "single_lassosum2" / "prs_method_final_res.txt",
        plot=PRS_METHOD_RUN_DIR / "single_lassosum2" / "prs_method_lassosum_plot.png",
        done=PRS_METHOD_RUN_DIR / "single_lassosum2.done",
    params:
        out_base=PRS_METHOD_RUN_DIR / "single_lassosum2",
        out_prefix=PRS_METHOD_RUN_DIR / "single_lassosum2" / "prs_method",
        ld_matrix_dir=lambda wildcards: prs_method_ld_matrix_dir("single_lassosum2"),
        afreq=PRS_AFREQ_FILE,
        script=PRS_SRC / "run_lassosum2.R",
    shell:
        """
        set -euo pipefail

        mkdir -p {params.out_base}
        cd {params.out_base}

        lassosum2_args="--rds {input.rds} --ss {input.ss} --bim {input.study_bim} --out {params.out_prefix} --ncores {threads} --pheno {input.pheno} --ld-matrix-dir {params.ld_matrix_dir}"
        if [[ -n "{params.afreq}" ]]; then
            lassosum2_args="$lassosum2_args --afreq {params.afreq}"
        fi

        Rscript {params.script} $lassosum2_args > {log} 2>&1

        touch {output.done}
        """


rule runSingleAncestryPRSice:
    log:
        OUT_DIR / "logs" / "runSingleAncestryPRSice.log",
    container:
        "oras://ghcr.io/mainsqu33ze/gdcgenomicsqc/prs:latest"
    conda:
        "../../envs/prs.yml"
    threads: 4
    resources:
        nodes=1,
        mem_mb=16000,
        runtime=240,
    input:
        resources=rules.preparePRSMethodResources.output.ready,
        ss=rules.alignSumstatsForPRS.output.aligned,
        pheno=rules.makeStudyPhenoFile.output.pheno,
        study_bed=rules.preparePRSInputs.output.study_bed,
        study_bim=rules.preparePRSInputs.output.study_bim,
        study_fam=rules.preparePRSInputs.output.study_fam,
    output:
        results=PRS_METHOD_RUN_DIR / "single_prsice" / "PRSice2" / "prs_method" / "PRSice2_outputs.prsice",
        best=PRS_METHOD_RUN_DIR / "single_prsice" / "PRSice2" / "prs_method" / "PRSice2_outputs.best",
        snps=PRS_METHOD_RUN_DIR / "single_prsice" / "PRSice2" / "prs_method" / "PRSice2_outputs.snps",
        done=PRS_METHOD_RUN_DIR / "single_prsice.done",
    params:
        study_prefix=PRS_STUDY_PLINK,
        out_base=PRS_METHOD_RUN_DIR / "single_prsice",
        out_prefix=PRS_METHOD_RUN_DIR / "single_prsice" / "PRSice2" / "prs_method",
        binary_target=PRS_CONFIG.get("binary_target", "F"),
        prsice_r=PRS_SRC / "PRSice.R",
        prsice_bin=PRS_SRC / "PRSice_linux",
        plink=PRS_PLINK_BIN,
        script=PRS_SRC / "run_PRSice2.sh",
    shell:
        """
        set -euo pipefail

        mkdir -p {params.out_base}
        export PRSICE_CMD="Rscript {params.prsice_r} --prsice {params.prsice_bin}"
        PLINK_BIN="{params.plink}"
        if [[ "$PLINK_BIN" != "plink" ]]; then
            export PLINK_CMD="$PLINK_BIN"
        fi

        bash {params.script} {input.ss} {params.study_prefix} {params.binary_target} {input.pheno} {params.out_base} {PRS_PIPELINE_PATH} {params.out_prefix} > {log} 2>&1

        touch {output.done}
        """


rule runSingleAncestryPRSCS:
    log:
        OUT_DIR / "logs" / "runSingleAncestryPRSCS.log",
    container:
        "oras://ghcr.io/mainsqu33ze/gdcgenomicsqc/prs:latest"
    conda:
        "../../envs/prs.yml"
    threads: 8
    resources:
        nodes=1,
        mem_mb=32000,
        runtime=720,
    input:
        resources=rules.preparePRSMethodResources.output.ready,
        ss=rules.alignSumstatsForPRS.output.aligned,
        study_bed=rules.preparePRSInputs.output.study_bed,
        study_bim=rules.preparePRSInputs.output.study_bim,
        study_fam=rules.preparePRSInputs.output.study_fam,
    output:
        combined=PRS_METHOD_RUN_DIR / "single_prscs" / "prs_pipeline" / "PRScs" / f"PRScs_{PRS_ANC1}_combined_weights.txt",
        score=PRS_METHOD_RUN_DIR / "single_prscs" / "prs_pipeline" / "PRScs" / f"PRScs_{PRS_ANC1}_score.sscore",
        rsq=PRS_METHOD_RUN_DIR / "single_prscs" / "prs_pipeline" / "PRScs" / f"{PRS_ANC1}_PRS_sscore_Rsqr.txt",
        adj_rsq=PRS_METHOD_RUN_DIR / "single_prscs" / "prs_pipeline" / "PRScs" / f"{PRS_ANC1}_adj_PRS_sscore_Rsqr.txt",
        done=PRS_METHOD_RUN_DIR / "single_prscs.done",
    params:
        study_prefix=PRS_STUDY_PLINK,
        out_base=PRS_METHOD_RUN_DIR / "single_prscs",
        anc=PRS_ANC1,
        plink2=PRS_PLINK2_BIN,
        path_code=PRS_PRSCS_PATH_CODE,
        ref_dir=PRS_PRSCS_REF_DIR,
        seed=PRS_PRSCS_SEED,
        path_python=lambda wildcards: PRS_PRSCS_PATH_PYTHON or "",
        rscript=lambda wildcards: "Rscript",
        script=PRS_SRC / "run_PRScs.sh",
    shell:
        """
        set -euo pipefail

        mkdir -p {params.out_base}
        CONFIG={params.out_base}/temp/PRScs_temp_config.txt
        mkdir -p "$(dirname "$CONFIG")"

        echo "target_sumstats_file={input.ss}" > "$CONFIG"
        echo "study_sample_plink={params.study_prefix}" >> "$CONFIG"
        echo "reference_SNPS_bim={params.study_prefix}" >> "$CONFIG"
        echo "output_dir={params.out_base}" >> "$CONFIG"
        echo "path_code={params.path_code}" >> "$CONFIG"
        echo "path_ref_dir={params.ref_dir}" >> "$CONFIG"
        echo "path_plink2={params.plink2}" >> "$CONFIG"
        echo "path_python={params.path_python}" >> "$CONFIG"
        echo "rscript={params.rscript}" >> "$CONFIG"
        echo "anc1={params.anc}" >> "$CONFIG"
        echo "seed={params.seed}" >> "$CONFIG"
        echo "prs_pipeline={PRS_PIPELINE_PATH}" >> "$CONFIG"

        bash {params.script} --c "$CONFIG" > {log} 2>&1

        touch {output.done}
        """


rule scoreTestPRS:
    log:
        OUT_DIR / "logs" / "scoreTestPRS.log",
    container:
        "oras://ghcr.io/mainsqu33ze/gdcgenomicsqc/prs:latest"
    conda:
        "../../envs/prs.yml"
    threads: 1
    resources:
        nodes=1,
        mem_mb=8000,
        runtime=720,
    input:
        resources=rules.preparePRSMethodResources.output.ready,
        ss=rules.alignSumstatsForPRS.output.aligned,
        ct_done=rules.runSingleAncestryCT.output.done,
        ct_results=rules.runSingleAncestryCT.output.results,
        ldpred2_done=rules.runSingleAncestryLDpred2.output.done,
        ldpred2_inf=rules.runSingleAncestryLDpred2.output.inf_weights,
        ldpred2_grid=rules.runSingleAncestryLDpred2.output.grid_weights,
        lassosum2_done=rules.runSingleAncestryLassosum2.output.done,
        lassosum2_weights=rules.runSingleAncestryLassosum2.output.weights,
        prsice2_done=rules.runSingleAncestryPRSice.output.done,
        prsice2_snps=rules.runSingleAncestryPRSice.output.snps,
        prscs_done=rules.runSingleAncestryPRSCS.output.done,
        prscs_combined=rules.runSingleAncestryPRSCS.output.combined,
        **PRS_TEST_INPUTS,
    output:
        ct_results=PRS_TEST_EVAL_DIR / "CT_results.txt",
        ct_scores=PRS_TEST_EVAL_DIR / "CT_scores.txt",
        ldpred2_inf_results=PRS_TEST_EVAL_DIR / "LDpred2_inf_results.txt",
        ldpred2_inf_scores=PRS_TEST_EVAL_DIR / "LDpred2_inf_scores.txt",
        ldpred2_grid_results=PRS_TEST_EVAL_DIR / "LDpred2_grid_results.txt",
        ldpred2_grid_scores=PRS_TEST_EVAL_DIR / "LDpred2_grid_scores.txt",
        lassosum2_results=PRS_TEST_EVAL_DIR / "lassosum2_results.txt",
        lassosum2_scores=PRS_TEST_EVAL_DIR / "lassosum2_scores.txt",
        prsice2_results=PRS_TEST_EVAL_DIR / "PRSice2_results.txt",
        prsice2_scores=PRS_TEST_EVAL_DIR / "PRSice2_scores.txt",
        prscsx_results=PRS_TEST_EVAL_DIR / f"PRScsx_{PRS_ANC1}_results.txt",
        prscsx_scores=PRS_TEST_EVAL_DIR / f"PRScsx_{PRS_ANC1}_scores.txt",
        done=PRS_METHOD_RUN_DIR / "scoreTestPRS.done",
    params:
        test_sample=PRS_TEST_SAMPLE,
        test_pca_args=(
            f"--test-pca-file {shlex.quote(str(PRS_TEST_PCA_FILE))}"
            if PRS_TEST_PCA_FILE
            else ""
        ),
        binary_target=PRS_CONFIG.get("binary_target", "F"),
        plink=PRS_PLINK_BIN,
        stage=PRS_TEST_STAGE_DIR,
        eval_dir=PRS_TEST_EVAL_DIR,
        ct_src=PRS_METHOD_RUN_DIR / "single_ct" / "CT",
        ldpred2_src=PRS_METHOD_RUN_DIR / "single_ldpred2",
        lassosum2_src=PRS_METHOD_RUN_DIR / "single_lassosum2",
        prsice2_src=PRS_METHOD_RUN_DIR / "single_prsice" / "PRSice2",
        prscs_label=PRS_ANC1,
        prsice_bin=PRS_SRC / "PRSice_linux",
        script=PRS_SRC / "score_test.sh",
    shell:
        """
        set -euo pipefail

        if [[ -z "{params.test_sample}" ]]; then
            echo "ERROR: prsPipeline.test_sample is not set; scoreTestPRS requires a test PLINK prefix." >&2
            exit 1
        fi
        if [[ ! -f "{params.test_sample}.fam" ]]; then
            echo "ERROR: test-sample PLINK files not found at {params.test_sample}.{{bed,bim,fam}}" >&2
            exit 1
        fi

        STAGE={params.stage}
        EVAL={params.eval_dir}
        mkdir -p "$STAGE" "$EVAL"

        # Bridge the native per-method layout to the upstream layout score_test.sh
        # expects ({{train_out_dir}}/CT, /LDpred2, /lassosum2, /PRSice2, /PRScsx).
        ln -sfn {params.ct_src}        "$STAGE/CT"
        ln -sfn {params.ldpred2_src}   "$STAGE/LDpred2"
        ln -sfn {params.lassosum2_src} "$STAGE/lassosum2"
        ln -sfn {params.prsice2_src}   "$STAGE/PRSice2"
        mkdir -p "$STAGE/PRScsx"
        ln -sfn {input.prscs_combined} "$STAGE/PRScsx/{params.prscs_label}_combined_weights.txt"

        # Provide PRSice and plink on PATH (score_test.sh calls the bare names).
        mkdir -p "$STAGE/bin"
        ln -sfn {params.prsice_bin} "$STAGE/bin/PRSice"
        export PATH="$STAGE/bin:$PATH"
        PLINK_BIN="{params.plink}"
        if [[ "$PLINK_BIN" != "plink" ]]; then
            export PATH="$(dirname "$PLINK_BIN"):$PATH"
        fi

        bash {params.script} \\
            --test-bfile {params.test_sample} \\
            --train-out-dir "$STAGE" \\
            {params.test_pca_args} \\
            --sumstats {input.ss} \\
            --path-repo {PRS_PIPELINE_PATH} \\
            --binary-flag {params.binary_target} \\
            --ran-ct true \\
            --ran-ldpred2 true \\
            --ran-lassosum2 true \\
            --ran-prsice2 true \\
            --ran-prscsx true > {log} 2>&1

        mv "$STAGE/test_evaluation"/* "$EVAL"/
        rm -rf "$STAGE"
        touch {output.done}
        """


CTSLEB_INPUT_DIR = Path(config.get("prsPipeline", {}).get("generated_input_dir", ""))
CTSLEB_PHENO_PREFIX = config.get("prsMethods", {}).get("multi_ctsleb", {}).get(
    "pheno_prefix",
    f"{PRS_ANC1}_HM3",
)

rule prepareCTSLEBSumstats:
    log:
        OUT_DIR / "logs" / "prepareCTSLEBSumstats.log",
    input:
        target_ss=rules.preparePRSInputs.output.target_sumstats,
        training_ss=rules.preparePRSInputs.output.training_sumstats,
        target_bim=rules.preparePRSInputs.output.study_bim,
        training_bim=rules.preparePRSInputs.output.study_anc2_bim,
    output:
        target_ss=CTSLEB_INPUT_DIR / "gwas" / "target_sumstats_CTSLEB.txt",
        training_ss=CTSLEB_INPUT_DIR / "gwas" / "training_sumstats_CTSLEB.txt",
    run:
        import pandas as pd

        def make_ctsleb_sumstats(ss_file, bim_file, out_file):
            ss = pd.read_csv(str(ss_file), sep="\t")
            bim = pd.read_csv(
                str(bim_file),
                sep=r"\s+",
                header=None,
                names=["CHR", "SNP", "CM", "BP", "A1_bim", "A2_bim"],
            )

            if "SNP" not in ss.columns:
                raise ValueError(f"Missing SNP column in {ss_file}")

            merged = ss.merge(bim[["SNP", "BP"]], on="SNP", how="left")
            missing_bp = merged["BP"].isna().sum()
            if missing_bp:
                raise ValueError(f"Missing BP for {missing_bp} rows in {out_file}")

            merged["BP"] = merged["BP"].astype(int)
            merged["rs_id"] = merged["SNP"]

            keep = ["CHR", "SNP", "BP", "A1", "A2", "BETA", "SE", "P", "N", "rs_id"]
            missing_cols = [col for col in keep if col not in merged.columns]
            if missing_cols:
                raise ValueError(
                    f"Missing columns in merged file {out_file}: {', '.join(missing_cols)}"
                )

            merged.loc[:, keep].to_csv(str(out_file), sep="\t", index=False)

        make_ctsleb_sumstats(input.target_ss, input.target_bim, output.target_ss)
        make_ctsleb_sumstats(input.training_ss, input.training_bim, output.training_ss)

        with open(str(log), "w", encoding="utf-8") as handle:
            handle.write(f"Wrote: {output.target_ss}\n")
            handle.write(f"Wrote: {output.training_ss}\n")

rule prepareCTSLEBPhenotypes:
    log:
        OUT_DIR / "logs" / "prepareCTSLEBPhenotypes.log",
    threads: 1
    resources:
        nodes=1,
        mem_mb=4000,
        runtime=60,
    input:
        fam=rules.preparePRSInputs.output.study_fam,
    output:
        tuning=CTSLEB_INPUT_DIR / "metadata" / f"{CTSLEB_PHENO_PREFIX}_tuning.pheno",
        validation=CTSLEB_INPUT_DIR / "metadata" / f"{CTSLEB_PHENO_PREFIX}_validation.pheno",
    shell:
        """
        set -euo pipefail

        mkdir -p $(dirname {output.tuning})

        N=$(wc -l < {input.fam})
        N_TUNE=$((N / 2))
        awk -v n_tune="$N_TUNE" -v tune="{output.tuning}" -v valid="{output.validation}" '
            NR <= n_tune {{ print $6 > tune; next }}
            {{ print $6 > valid }}
        ' {input.fam}

        echo "Created CTSLEB phenotype splits:" > {log}
        wc -l {output.tuning} {output.validation} >> {log}
        """

rule runMultiAncestryPRSCSx:
    log:
        OUT_DIR / "logs" / "runMultiAncestryPRSCSx.log",
    container:
        "oras://ghcr.io/mainsqu33ze/gdcgenomicsqc/prs:latest"
    conda:
        "../../envs/prs.yml"
    threads: 4
    resources:
        nodes=1,
        mem_mb=32000,
        runtime=720,
    input:
        resources=rules.preparePRSMethodResources.output.ready,
        ss_target=rules.preparePRSInputs.output.target_sumstats,
        ss_training=rules.preparePRSInputs.output.training_sumstats,
        study_bed=rules.preparePRSInputs.output.study_bed,
        study_bim=rules.preparePRSInputs.output.study_bim,
        study_fam=rules.preparePRSInputs.output.study_fam,
        study_anc2_bed=rules.preparePRSInputs.output.study_anc2_bed,
        study_anc2_bim=rules.preparePRSInputs.output.study_anc2_bim,
        study_anc2_fam=rules.preparePRSInputs.output.study_anc2_fam,
    output:
        combined_anc1=PRS_METHOD_RUN_DIR / "multi_prscsx" / "prs_pipeline" / "PRScsx" / f"PRScsx_{PRS_ANC1}_combined_weights.txt",
        combined_anc2=PRS_METHOD_RUN_DIR / "multi_prscsx" / "prs_pipeline" / "PRScsx" / f"PRScsx_{PRS_ANC2}_combined_weights.txt",
        score_anc1=PRS_METHOD_RUN_DIR / "multi_prscsx" / "prs_pipeline" / "PRScsx" / f"PRScsx_joint_{PRS_ANC1}_score.sscore",
        score_anc2=PRS_METHOD_RUN_DIR / "multi_prscsx" / "prs_pipeline" / "PRScsx" / f"PRScsx_joint_{PRS_ANC2}_score.sscore",
        rsq_anc1=PRS_METHOD_RUN_DIR / "multi_prscsx" / "prs_pipeline" / "PRScsx" / f"{PRS_ANC1}_PRS_sscore_Rsqr.txt",
        adj_rsq_anc1=PRS_METHOD_RUN_DIR / "multi_prscsx" / "prs_pipeline" / "PRScsx" / f"{PRS_ANC1}_adj_PRS_sscore_Rsqr.txt",
        rsq_anc2=PRS_METHOD_RUN_DIR / "multi_prscsx" / "prs_pipeline" / "PRScsx" / f"{PRS_ANC2}_PRS_sscore_Rsqr.txt",
        adj_rsq_anc2=PRS_METHOD_RUN_DIR / "multi_prscsx" / "prs_pipeline" / "PRScsx" / f"{PRS_ANC2}_adj_PRS_sscore_Rsqr.txt",
        done=PRS_METHOD_RUN_DIR / "multi_prscsx.done",
    params:
        study_prefix=PRS_STUDY_PLINK,
        study_anc2_prefix=PRS_OUT_DIR / "anc2_plink_files" / f"{PRS_ANC2}_simulation_study_sample",
        out_base=PRS_METHOD_RUN_DIR / "multi_prscsx",
        anc1=PRS_ANC1,
        anc2=PRS_ANC2,
        plink2=PRS_PLINK2_BIN,
        path_code=PRS_PRSCSX_PATH_CODE,
        ref_dir=PRS_PRSCSX_REF_DIR,
        seed=PRS_PRSCSX_SEED,
        path_python=lambda wildcards: PRS_PRSCSX_PATH_PYTHON or "",
        rscript=lambda wildcards: "Rscript",
        script=PRS_SRC / "run_PRScsx.sh",
    shell:
        """
        set -euo pipefail

        mkdir -p {params.out_base}
        CONFIG={params.out_base}/temp/PRScsx_temp_config.txt
        mkdir -p "$(dirname "$CONFIG")"

        echo "path_data_root={params.out_base}" > "$CONFIG"
        echo "target_sumstats_file={input.ss_target}" >> "$CONFIG"
        echo "training_sumstats_file={input.ss_training}" >> "$CONFIG"
        echo "study_sample_plink={params.study_prefix}" >> "$CONFIG"
        echo "study_sample_plink_anc2={params.study_anc2_prefix}" >> "$CONFIG"
        echo "reference_SNPS_bim={params.study_prefix}" >> "$CONFIG"
        echo "output_dir={params.out_base}" >> "$CONFIG"
        echo "path_code={params.path_code}" >> "$CONFIG"
        echo "path_ref_dir={params.ref_dir}" >> "$CONFIG"
        echo "path_plink2={params.plink2}" >> "$CONFIG"
        echo "path_python={params.path_python}" >> "$CONFIG"
        echo "rscript={params.rscript}" >> "$CONFIG"
        echo "anc1={params.anc1}" >> "$CONFIG"
        echo "anc2={params.anc2}" >> "$CONFIG"
        echo "seed={params.seed}" >> "$CONFIG"
        echo "prs_pipeline={PRS_PIPELINE_PATH}" >> "$CONFIG"

        bash {params.script} --c "$CONFIG" > {log} 2>&1

        touch {output.done}
        """


rule runMultiAncestryLDpred2:
    log:
        OUT_DIR / "logs" / "runMultiAncestryLDpred2.log",
    threads: 4
    resources:
        nodes=1,
        mem_mb=64000,
        runtime=720,
    input:
        resources=rules.preparePRSMethodResources.output.ready,
        env=rules.preparePRSInputs.output.env,
        target_sumstats=rules.preparePRSInputs.output.target_sumstats,
        training_sumstats=rules.preparePRSInputs.output.training_sumstats,
        study_bed=rules.preparePRSInputs.output.study_bed,
        study_anc2_bed=rules.preparePRSInputs.output.study_anc2_bed,
        study_pheno=rules.preparePRSInputs.output.target_study_pheno,
        study_anc2_pheno=rules.preparePRSInputs.output.training_study_pheno,
    output:
        done=PRS_METHOD_RUN_DIR / "multi_ldpred2.done",
    params:
        method="multi_ldpred2",
        command=prs_method_command_quoted("multi_ldpred2"),
        extra=prs_method_extra_args("multi_ldpred2"),
        out_dir=PRS_METHOD_RUN_DIR / "multi_ldpred2",
        script=Path(workflow.basedir) / "scripts" / "run_prs_pipeline_adapter.sh",
    shell:
        """
        PRS_METHOD_COMMAND={params.command} bash {params.script} \
            --method {params.method} \
            --prs-inputs-env {input.env} \
            --resource-dir {PRS_RESOURCE_DIR} \
            --out-dir {params.out_dir} \
            {params.extra} \
            --done {output.done} \
            > {log} 2>&1
        """



rule runMultiAncestrySDPRS:
    log:
        OUT_DIR / "logs" / "runMultiAncestrySDPRS.log",
    threads: 4
    resources:
        nodes=1,
        mem_mb=64000,
        runtime=720,
    input:
        resources=rules.preparePRSMethodResources.output.ready,
        env=rules.preparePRSInputs.output.env,
        target_sumstats=rules.preparePRSInputs.output.target_sumstats,
        training_sumstats=rules.preparePRSInputs.output.training_sumstats,
        study_pheno=rules.preparePRSInputs.output.target_study_pheno,
        study_anc2_pheno=rules.preparePRSInputs.output.training_study_pheno,
    output:
        done=PRS_METHOD_RUN_DIR / "multi_sdprs.done",
    params:
        method="multi_sdprs",
        command=prs_method_command_quoted("multi_sdprs"),
        extra=prs_method_extra_args("multi_sdprs"),
        out_dir=PRS_METHOD_RUN_DIR / "multi_sdprs",
        script=Path(workflow.basedir) / "scripts" / "run_prs_pipeline_adapter.sh",
    shell:
        """
        PRS_METHOD_COMMAND={params.command} bash {params.script} \
            --method {params.method} \
            --prs-inputs-env {input.env} \
            --resource-dir {PRS_RESOURCE_DIR} \
            --out-dir {params.out_dir} \
            {params.extra} \
            --done {output.done} \
            > {log} 2>&1
        """
