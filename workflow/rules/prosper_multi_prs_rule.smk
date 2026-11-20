PROSPER_INPUT_DIR = Path(config.get("prsPipeline", {}).get("generated_input_dir", ""))

# PROSPER (joint multi-ancestry PRS: penalized regression + super-learning) runs
# inside the upstream `singleprshelper` container (oras://ghcr.io/mainsqu33ze/
# gdcgenomicsqc/singleprshelper:latest) — PROSPER.R/lassosum2.R need
# RcppArmadillo, doMC, etc. + a compiler (sourceCpp) that the generic `prs`
# image does not ship. The PROSPER R package lives on the host/gdc (bound via
# /projects in the profiles); a writable per-run copy carrying ref_bim.txt is
# made by preparePROSPERPackage because PROSPER.R/lassosum2.R read
# <PATH_package>/ref_bim.txt and the shared checkout is non-writable.
PROSPER_METHOD_CONFIG = config.get("prsMethods", {}).get("multi_prosper", {})
PROSPER_SOFTWARE_DIR = PROSPER_METHOD_CONFIG.get(
    "software_dir",
    "/projects/standard/gdc/public/prs_methods/scripts/PROSPER",
)
PROSPER_TARGET_POP = PROSPER_METHOD_CONFIG.get("target_pop", PRS_ANC1)
PROSPER_AUX_POP = PROSPER_METHOD_CONFIG.get("aux_pop", PRS_ANC2)
PROSPER_CHROM = PROSPER_METHOD_CONFIG.get("chrom", "1-22")
PROSPER_NCORES = PROSPER_METHOD_CONFIG.get("ncores", 1)
PROSPER_LL = PROSPER_METHOD_CONFIG.get("ll", 5)
PROSPER_LC = PROSPER_METHOD_CONFIG.get("lc", 5)
PROSPER_CONTAINER = "oras://ghcr.io/mainsqu33ze/gdcgenomicsqc/singleprshelper:latest"
# Writable per-run package copy + ref_bim.txt (see preparePROSPERPackage).
PROSPER_PKG_DIR = PRS_METHOD_RUN_DIR / "multi_prosper" / "PROSPER_pkg"
# Where PROSPER.R / tuning_testing.R write (after_ensemble_<TARGET_POP>/…).
PROSPER_RUN_DIR = PRS_METHOD_RUN_DIR / "multi_prosper" / "PROSPER"
PROSPER_AFTER_ENSEMBLE_DIR = PROSPER_RUN_DIR / ("after_ensemble_" + PROSPER_TARGET_POP)
PROSPER_PLINK2 = (
    PROSPER_METHOD_CONFIG.get("plink2")
    or config.get("prsMethods", {}).get("plink2")
    or "plink2"
)
PROSPER_RSCRIPT = PROSPER_METHOD_CONFIG.get("rscript", "Rscript")



rule preparePROSPERSumstats:
    log:
        OUT_DIR / "logs" / "preparePROSPERSumstats.log",
    input:
        target_ss=rules.preparePRSInputs.output.target_sumstats,
        training_ss=rules.preparePRSInputs.output.training_sumstats,
        target_bim=rules.preparePRSInputs.output.study_bim,
    output:
        target_ss=PROSPER_INPUT_DIR / "gwas" / "target_sumstats_PROSPER.txt",
        training_ss=PROSPER_INPUT_DIR / "gwas" / "training_sumstats_PROSPER.txt",
    run:
        import pandas as pd

        def make_prosper(ss_file, bim_file, out_file):
            ss = pd.read_csv(str(ss_file), sep="\t")
            ss.columns = [c.lower() for c in ss.columns]

            rename = {
                "snp": "rsid",
                "a2": "a0",
                "se": "beta_se",
                "n": "n_eff",
                "bp": "pos",
            }
            ss = ss.rename(columns={k: v for k, v in rename.items() if k in ss.columns})

            bim = pd.read_csv(
                str(bim_file),
                sep=r"\s+",
                header=None,
                names=["chr_bim", "rsid", "cm", "pos", "a1_bim", "a0_bim"],
            )

            if "pos" not in ss.columns:
                ss = ss.merge(bim[["rsid", "pos"]], on="rsid", how="left")

            keep = ["rsid", "chr", "pos", "a1", "a0", "beta", "beta_se", "p", "n_eff"]
            missing = [c for c in keep if c not in ss.columns]
            if missing:
                raise ValueError(f"Missing columns in {ss_file}: {', '.join(missing)}")

            ss = ss.loc[ss["rsid"].notna() & ss["beta"].notna() & ss["beta_se"].notna(), keep]
            ss.to_csv(str(out_file), sep="\t", index=False)

        make_prosper(input.target_ss, input.target_bim, output.target_ss)
        make_prosper(input.training_ss, input.target_bim, output.training_ss)

        with open(str(log), "w", encoding="utf-8") as handle:
            handle.write(f"Wrote: {output.target_ss}\n")
            handle.write(f"Wrote: {output.training_ss}\n")


# PROSPER.R / lassosum2.R hard-code <PATH_package>/ref_bim.txt + scripts/ and
# <PATH_package>/<pop>/chr*_LD.RData. The package's LD reference dirs (AFR,
# EUR, ...) are large (multi-GB each), so the per-run package dir mirrors the
# shared checkout with symlinks for those, carrying a real ref_bim.txt
# (preferring the package's own — its LD blocks were built against it).
rule preparePROSPERPackage:
    log:
        OUT_DIR / "logs" / "preparePROSPERPackage.log",
    threads: 1
    resources:
        nodes=1,
        mem_mb=2000,
        runtime=60,
    input:
        bim=rules.preparePRSInputs.output.study_bim,
    output:
        ref_bim=PROSPER_PKG_DIR / "ref_bim.txt",
        prosper_r=PROSPER_PKG_DIR / "scripts" / "PROSPER.R",
        lassosum2_r=PROSPER_PKG_DIR / "scripts" / "lassosum2.R",
        tuning_testing_r=PROSPER_PKG_DIR / "scripts" / "tuning_testing.R",
    params:
        software_dir=PROSPER_SOFTWARE_DIR,
    shell:
        """
        set -euo pipefail

        if [[ ! -d {params.software_dir}/scripts ]]; then
            echo "PROSPER package not found at {params.software_dir}/scripts" >&2
            echo "Set prsMethods.multi_prosper.software_dir to the PROSPER R package" >&2
            exit 1
        fi

        PKG_DIR=$(dirname {output.ref_bim})
        mkdir -p "$PKG_DIR/scripts"
        cp -r {params.software_dir}/scripts/. "$PKG_DIR/scripts/"

        # plink2 --score refuses to impute allele frequencies when the scoring
        # set has <50 samples. Patch the vendored copies (not the shared
        # checkout) so every --score call loads --read-freq <bfile>.afreq,
        # which preparePROSPERTuningTestingSets generates alongside each set.
        sed -i \
            -e 's|" --bfile ", bfile_tuning,|" --bfile ", bfile_tuning,\\n                " --read-freq ", paste0(bfile_tuning, ".afreq"),|' \
            -e 's|" --bfile ", bfile_testing,|" --bfile ", bfile_testing,\\n                " --read-freq ", paste0(bfile_testing, ".afreq"),|' \
            "$PKG_DIR/scripts/lassosum2.R"
        sed -i \
            -e 's|" --bfile ",opt$bfile_tuning,|" --bfile ",opt$bfile_tuning,\\n              " --read-freq ", paste0(opt$bfile_tuning, ".afreq"),|' \
            -e 's|" --bfile ",opt$bfile_testing,|" --bfile ",opt$bfile_testing,\\n              " --read-freq ", paste0(opt$bfile_testing, ".afreq"),|' \
            "$PKG_DIR/scripts/tuning_testing.R"

        # plink2 --score emits 5 baseline/aux .sscore columns
        # (#FID IID PHENO1 ALLELE_CT NAMED_ALLELE_DOSAGE_SUM) before the score
        # _SUM columns; the upstream scripts drop only 4, leaving one extra
        # soup column that makes the ensemble matmult non-conformable. Drop 5.
        sed -i 's|SCORE <- SCORE\[,-1:-4\]|SCORE <- SCORE[,-1:-5]|g' \
            "$PKG_DIR/scripts/lassosum2.R" \
            "$PKG_DIR/scripts/tuning_testing.R"

        if [[ -f {params.software_dir}/ref_bim.txt ]]; then
            # Reference BIM that the shipped LD blocks were built against.
            cp {params.software_dir}/ref_bim.txt {output.ref_bim}
        else
            # Fallback: derive a 6-col BIM from the study (no matching LD blocks).
            awk 'BEGIN{{OFS="\\t"}} {{print $1, $2, $3, $4, $5, $6}}' {input.bim} > {output.ref_bim}
        fi

        # Symlink per-population LD reference dirs (read-only loads inside the
        # container); skip non-LD artifacts shipped in the package root.
        for pop_dir in {params.software_dir}/*/; do
            pop=$(basename "$pop_dir")
            case "$pop" in
                scripts|.git|example|PROSPER_example_results|PROSPER_example_results_2) continue ;;
            esac
            if compgen -G "$pop_dir"chr*_LD.RData > /dev/null; then
                ln -sfn "$pop_dir" "$PKG_DIR/$pop"
            fi
        done

        echo "Copied PROSPER package from {params.software_dir}" > {log}
        echo "Package dir: $PKG_DIR" >> {log}
        echo "Wrote {output.ref_bim} ($(wc -l < {output.ref_bim}) SNPs)" >> {log}
        echo "LD ref dirs:" >> {log}
        find "$PKG_DIR" -maxdepth 1 -type l | sort >> {log}
        echo "LD block files: $(find -L "$PKG_DIR" -name 'chr*_LD.RData' | wc -l)" >> {log}
        """



rule preparePROSPERTuningTestingSets:
    log:
        OUT_DIR / "logs" / "preparePROSPERTuningTestingSets.log",
    threads: 1
    resources:
        nodes=1,
        mem_mb=8000,
        runtime=120,
    input:
        bed=rules.preparePRSInputs.output.study_bed,
        bim=rules.preparePRSInputs.output.study_bim,
        fam=rules.preparePRSInputs.output.study_fam,
        aux_bed=rules.preparePRSInputs.output.study_anc2_bed,
        aux_fam=rules.preparePRSInputs.output.study_anc2_fam,
    output:
        tuning_keep=PROSPER_INPUT_DIR / "metadata" / f"{PRS_ANC1}_PROSPER_tuning.keep",
        testing_keep=PROSPER_INPUT_DIR / "metadata" / f"{PRS_ANC1}_PROSPER_testing.keep",
        tuning_pheno=PROSPER_INPUT_DIR / "metadata" / f"{PRS_ANC1}_PROSPER_tuning.pheno",
        testing_pheno=PROSPER_INPUT_DIR / "metadata" / f"{PRS_ANC1}_PROSPER_testing.pheno",
        tuning_bed=PROSPER_INPUT_DIR / "prosper" / f"{PRS_ANC1}_tuning.bed",
        tuning_bim=PROSPER_INPUT_DIR / "prosper" / f"{PRS_ANC1}_tuning.bim",
        tuning_fam=PROSPER_INPUT_DIR / "prosper" / f"{PRS_ANC1}_tuning.fam",
        testing_bed=PROSPER_INPUT_DIR / "prosper" / f"{PRS_ANC1}_testing.bed",
        testing_bim=PROSPER_INPUT_DIR / "prosper" / f"{PRS_ANC1}_testing.bim",
        testing_fam=PROSPER_INPUT_DIR / "prosper" / f"{PRS_ANC1}_testing.fam",
        aux_tuning_bed=PROSPER_INPUT_DIR / "prosper" / f"{PRS_ANC2}_tuning.bed",
        aux_tuning_bim=PROSPER_INPUT_DIR / "prosper" / f"{PRS_ANC2}_tuning.bim",
        aux_tuning_fam=PROSPER_INPUT_DIR / "prosper" / f"{PRS_ANC2}_tuning.fam",
        aux_tuning_pheno=PROSPER_INPUT_DIR / "metadata" / f"{PRS_ANC2}_PROSPER_tuning.pheno",
        # plink2 --score refuses to impute allele frequencies from <50 samples;
        # the vendored lassosum2.R / tuning_testing.R are patched (see
        # preparePROSPERPackage) to pass --read-freq pointing at these.
        tuning_afreq=PROSPER_INPUT_DIR / "prosper" / f"{PRS_ANC1}_tuning.afreq",
        testing_afreq=PROSPER_INPUT_DIR / "prosper" / f"{PRS_ANC1}_testing.afreq",
        aux_tuning_afreq=PROSPER_INPUT_DIR / "prosper" / f"{PRS_ANC2}_tuning.afreq",
    params:
        plink2=config.get("prsMethods", {}).get("plink2", "/projects/standard/gdc/public/plink2"),
        study_prefix=lambda wc, input: str(input.bed)[:-4],
        tuning_prefix=lambda wc, output: str(output.tuning_bed)[:-4],
        testing_prefix=lambda wc, output: str(output.testing_bed)[:-4],
        aux_prefix=lambda wc, output: str(output.aux_tuning_bed)[:-4],
        aux_base=lambda wc, input: str(input.aux_bed)[:-4],
    shell:
        """
        set -euo pipefail

        mkdir -p $(dirname {output.tuning_keep})
        mkdir -p $(dirname {output.tuning_bed})

        N=$(wc -l < {input.fam})
        N_TUNE=$((N / 2))
        N_TEST=$((N - N_TUNE))

        awk '{{print $1, $2}}' {input.fam} | head -n $N_TUNE > {output.tuning_keep}
        awk '{{print $1, $2}}' {input.fam} | tail -n $N_TEST > {output.testing_keep}
        awk '{{print $1, $2, $6}}' {input.fam} | head -n $N_TUNE > {output.tuning_pheno}
        awk '{{print $1, $2, $6}}' {input.fam} | tail -n $N_TEST > {output.testing_pheno}

        {params.plink2} \
          --bfile {params.study_prefix} \
          --keep {output.tuning_keep} \
          --make-bed \
          --freq \
          --out {params.tuning_prefix}

        {params.plink2} \
          --bfile {params.study_prefix} \
          --keep {output.testing_keep} \
          --make-bed \
          --freq \
          --out {params.testing_prefix}

        # Aux (ANC2) tuning set for lassosum2: whole ANC2 study sample + fam pheno.
        {params.plink2} \
          --bfile {params.aux_base} \
          --make-bed \
          --freq \
          --out {params.aux_prefix}
        awk '{{print $1, $2, $6}}' {input.aux_fam} > {output.aux_tuning_pheno}

        echo "Prepared PROSPER tuning/testing PLINK sets" > {log}
        wc -l {output.tuning_keep} {output.testing_keep} >> {log}
        wc -l {output.tuning_pheno} {output.testing_pheno} >> {log}
        wc -l {input.aux_fam} >> {log}
        """


rule preparePROSPERLassosumParams:
    log:
        OUT_DIR / "logs" / "preparePROSPERLassosumParams.log",
    threads: 4
    resources:
        nodes=1,
        mem_mb=32000,
        runtime=720,
    container:
        PROSPER_CONTAINER
    input:
        ref_bim=rules.preparePROSPERPackage.output.ref_bim,
        target_sumstats=rules.preparePROSPERSumstats.output.target_ss,
        training_sumstats=rules.preparePROSPERSumstats.output.training_ss,
        tuning_bed=rules.preparePROSPERTuningTestingSets.output.tuning_bed,
        tuning_bim=rules.preparePROSPERTuningTestingSets.output.tuning_bim,
        tuning_fam=rules.preparePROSPERTuningTestingSets.output.tuning_fam,
        tuning_pheno=rules.preparePROSPERTuningTestingSets.output.tuning_pheno,
        tuning_afreq=rules.preparePROSPERTuningTestingSets.output.tuning_afreq,
        aux_tuning_bed=rules.preparePROSPERTuningTestingSets.output.aux_tuning_bed,
        aux_tuning_pheno=rules.preparePROSPERTuningTestingSets.output.aux_tuning_pheno,
        aux_tuning_afreq=rules.preparePROSPERTuningTestingSets.output.aux_tuning_afreq,
    output:
        target_param=PRS_METHOD_RUN_DIR / "multi_prosper" / "lassosum2_AFR" / "param.txt",
        aux_param=PRS_METHOD_RUN_DIR / "multi_prosper" / "lassosum2_EUR" / "param.txt",
    params:
        package=PROSPER_PKG_DIR,
        plink2=PROSPER_PLINK2,
        rscript=PROSPER_RSCRIPT,
        target_pop=PROSPER_TARGET_POP,
        aux_pop=PROSPER_AUX_POP,
        chrom=PROSPER_CHROM,
        ncores=PROSPER_NCORES,
        tuning_prefix=lambda wc, input: str(input.tuning_bed)[:-4],
        aux_prefix=lambda wc, input: str(input.aux_tuning_bed)[:-4],
        target_out=lambda wc, output: str(Path(output.target_param).parent),
        aux_out=lambda wc, output: str(Path(output.aux_param).parent),
    shell:
        r"""
        set -euo pipefail

        # NOTE: no CC/CXX/MAKEFLAGS overrides here. singleprshelper is a conda
        # R image whose R sees itself compiled with the conda toolchain
        # (x86_64-conda-linux-gnu-cc/c++), and Rcpp::sourceCpp compiles cleanly
        # with R's own Makeconf. Pointing CC at a bare "gcc" would break it.
        mkdir -p "$(dirname {output.target_param})" "$(dirname {output.aux_param})"

        find_lassosum_param() {{
            local search_dir="$1"
            local dest="$2"
            local matches_file
            local n_found
            local found

            matches_file="${{search_dir}}/.lassosum_param_candidates.txt"
            find "$search_dir" -type f \
                ! -path '*/tmp/*' \
                ! -name '.lassosum_param_candidates.txt' \
                -print0 \
              | xargs -0 grep -l -E '(^|[[:space:]])delta0([[:space:]]|$).*lambda0|(^|[[:space:]])lambda0([[:space:]]|$).*delta0' \
              > "$matches_file" 2>/dev/null || true

            n_found=$(wc -l < "$matches_file")
            if [[ "$n_found" -eq 0 ]]; then
                echo "Could not find lassosum2 parameter file with delta0/lambda0 under $search_dir" >&2
                find "$search_dir" -maxdepth 3 -type f | sort >&2
                exit 1
            fi
            if [[ "$n_found" -gt 1 ]]; then
                echo "Found multiple lassosum2 parameter candidates under $search_dir; refusing to guess." >&2
                cat "$matches_file" >&2
                exit 1
            fi

            found=$(cat "$matches_file")
            cp "$found" "$dest"
            echo "Copied lassosum2 parameter file: $found -> $dest"
        }}

        echo "Running PROSPER simplified lassosum2 for {params.target_pop}" > {log}
        {params.rscript} {params.package}/scripts/lassosum2.R \
          --PATH_package {params.package} \
          --PATH_out {params.target_out} \
          --PATH_plink {params.plink2} \
          --FILE_sst {input.target_sumstats} \
          --pop {params.target_pop} \
          --chrom {params.chrom} \
          --bfile_tuning {params.tuning_prefix} \
          --pheno_tuning {input.tuning_pheno} \
          --NCORES {params.ncores} \
          >> {log} 2>&1
        find_lassosum_param "{params.target_out}" "{output.target_param}" >> {log} 2>&1

        echo "Running PROSPER simplified lassosum2 for {params.aux_pop}" >> {log}
        {params.rscript} {params.package}/scripts/lassosum2.R \
          --PATH_package {params.package} \
          --PATH_out {params.aux_out} \
          --PATH_plink {params.plink2} \
          --FILE_sst {input.training_sumstats} \
          --pop {params.aux_pop} \
          --chrom {params.chrom} \
          --bfile_tuning {params.aux_prefix} \
          --pheno_tuning {input.aux_tuning_pheno} \
          --NCORES {params.ncores} \
          >> {log} 2>&1
        find_lassosum_param "{params.aux_out}" "{output.aux_param}" >> {log} 2>&1
        """


rule runMultiAncestryPROSPER:
    log:
        OUT_DIR / "logs" / "runMultiAncestryPROSPER.log",
    threads: 4
    resources:
        nodes=1,
        mem_mb=32000,
        runtime=720,
    container:
        PROSPER_CONTAINER
    input:
        ref_bim=rules.preparePROSPERPackage.output.ref_bim,
        target_sumstats=rules.preparePROSPERSumstats.output.target_ss,
        training_sumstats=rules.preparePROSPERSumstats.output.training_ss,
        lassosum_target=rules.preparePROSPERLassosumParams.output.target_param,
        lassosum_aux=rules.preparePROSPERLassosumParams.output.aux_param,
        tuning_bed=rules.preparePROSPERTuningTestingSets.output.tuning_bed,
        tuning_pheno=rules.preparePROSPERTuningTestingSets.output.tuning_pheno,
        tuning_afreq=rules.preparePROSPERTuningTestingSets.output.tuning_afreq,
        testing_bed=rules.preparePROSPERTuningTestingSets.output.testing_bed,
        testing_pheno=rules.preparePROSPERTuningTestingSets.output.testing_pheno,
        testing_afreq=rules.preparePROSPERTuningTestingSets.output.testing_afreq,
        wrapper=Path(workflow.basedir) / "scripts" / "run_prosper_wrapper.R",
    output:
        done=PRS_METHOD_RUN_DIR / "multi_prosper.done",
        r2=PROSPER_AFTER_ENSEMBLE_DIR / "R2.txt",
        prs=PROSPER_AFTER_ENSEMBLE_DIR / "PROSPER_prs_file.txt",
    params:
        package=PROSPER_PKG_DIR,
        plink2=PROSPER_PLINK2,
        rscript=PROSPER_RSCRIPT,
        target_pop=PROSPER_TARGET_POP,
        aux_pop=PROSPER_AUX_POP,
        chrom=PROSPER_CHROM,
        ncores=PROSPER_NCORES,
        ll=PROSPER_LL,
        lc=PROSPER_LC,
        out_dir=PROSPER_RUN_DIR,
        tuning_prefix=lambda wc, input: str(input.tuning_bed)[:-4],
        testing_prefix=lambda wc, input: str(input.testing_bed)[:-4],
    shell:
        """
        set -euo pipefail
        mkdir -p {params.out_dir}
        {params.rscript} {input.wrapper} \\
            --package {params.package} \\
            --out-dir {params.out_dir} \\
            --plink2 {params.plink2} \\
            --target-sumstats {input.target_sumstats} \\
            --aux-sumstats {input.training_sumstats} \\
            --pop {params.target_pop},{params.aux_pop} \\
            --lassosum-param-files {input.lassosum_target},{input.lassosum_aux} \\
            --prefix {params.target_pop} \\
            --tuning-bfile {params.tuning_prefix} \\
            --testing-bfile {params.testing_prefix} \\
            --pheno-tuning {input.tuning_pheno} \\
            --pheno-testing {input.testing_pheno} \\
            --chrom {params.chrom} \\
            --ll {params.ll} --lc {params.lc} --ncores {params.ncores} \\
            > {log} 2>&1
        touch {output.done}
        """