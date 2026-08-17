PROSPER_INPUT_DIR = Path(config.get("prsPipeline", {}).get("generated_input_dir", ""))






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
    output:
        tuning_keep=PROSPER_INPUT_DIR / "metadata" / "AFR_PROSPER_tuning.keep",
        testing_keep=PROSPER_INPUT_DIR / "metadata" / "AFR_PROSPER_testing.keep",
        tuning_pheno=PROSPER_INPUT_DIR / "metadata" / "AFR_PROSPER_tuning.pheno",
        testing_pheno=PROSPER_INPUT_DIR / "metadata" / "AFR_PROSPER_testing.pheno",
        tuning_bed=PROSPER_INPUT_DIR / "prosper" / "AFR_tuning.bed",
        tuning_bim=PROSPER_INPUT_DIR / "prosper" / "AFR_tuning.bim",
        tuning_fam=PROSPER_INPUT_DIR / "prosper" / "AFR_tuning.fam",
        testing_bed=PROSPER_INPUT_DIR / "prosper" / "AFR_testing.bed",
        testing_bim=PROSPER_INPUT_DIR / "prosper" / "AFR_testing.bim",
        testing_fam=PROSPER_INPUT_DIR / "prosper" / "AFR_testing.fam",
        ref_bim="/scratch.global/saonli/GDCGenomicsQC/prs_resources/software/PROSPER/ref_bim.txt",
    params:
        plink2=config.get("prsMethods", {}).get("plink2", "/projects/standard/gdc/public/plink2"),
        study_prefix=lambda wc, input: str(input.bed)[:-4],
        tuning_prefix=lambda wc, output: str(output.tuning_bed)[:-4],
        testing_prefix=lambda wc, output: str(output.testing_bed)[:-4],
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
          --out {params.tuning_prefix}

        {params.plink2} \
          --bfile {params.study_prefix} \
          --keep {output.testing_keep} \
          --make-bed \
          --out {params.testing_prefix}

        awk 'BEGIN{{OFS="\t"}} {{print $1,$2,$3,$4,$5,$6}}' \
        {output.tuning_bim} > {output.ref_bim}
        echo "Prepared PROSPER tuning/testing PLINK sets" > {log}
        wc -l {output.tuning_keep} {output.testing_keep} >> {log}
        wc -l {output.tuning_pheno} {output.testing_pheno} >> {log}
        """


rule preparePROSPERLassosumParams:
    log:
        OUT_DIR / "logs" / "preparePROSPERLassosumParams.log",
    threads: 4
    resources:
        nodes=1,
        mem_mb=32000,
        runtime=720,
    input:
        target_sumstats=rules.preparePROSPERSumstats.output.target_ss,
        training_sumstats=rules.preparePROSPERSumstats.output.training_ss,
        tuning_bed=rules.preparePROSPERTuningTestingSets.output.tuning_bed,
        tuning_bim=rules.preparePROSPERTuningTestingSets.output.tuning_bim,
        tuning_fam=rules.preparePROSPERTuningTestingSets.output.tuning_fam,
        tuning_pheno=rules.preparePROSPERTuningTestingSets.output.tuning_pheno,
    output:
        target_param=PRS_METHOD_RUN_DIR / "multi_prosper" / "lassosum2_AFR" / "param.txt",
        aux_param=PRS_METHOD_RUN_DIR / "multi_prosper" / "lassosum2_EUR" / "param.txt",
    params:
        package=config.get("prsMethods", {}).get("multi_prosper", {}).get(
            "software_dir",
            "/scratch.global/saonli/GDCGenomicsQC/prs_resources/software/PROSPER",
        ),
        plink2=config.get("prsMethods", {}).get("multi_prosper", {}).get(
            "plink2",
            config.get("prsMethods", {}).get("plink2", "/projects/standard/gdc/public/plink2"),
        ),
        rscript=config.get("prsMethods", {}).get("multi_prosper", {}).get(
            "rscript",
            "/scratch.global/saonli/conda-envs/ctsleb/bin/Rscript",
        ),
        target_pop=config.get("prsMethods", {}).get("multi_prosper", {}).get("target_pop", "AFR"),
        aux_pop=config.get("prsMethods", {}).get("multi_prosper", {}).get("aux_pop", "EUR"),
        chrom=config.get("prsMethods", {}).get("multi_prosper", {}).get("chrom", "1-22"),
        ncores=config.get("prsMethods", {}).get("multi_prosper", {}).get("ncores", 1),
        tuning_prefix=lambda wc, input: str(input.tuning_bed)[:-4],
        target_out=lambda wc, output: str(Path(output.target_param).parent),
        aux_out=lambda wc, output: str(Path(output.aux_param).parent),
    shell:
        r"""
        set -euo pipefail
        module load gcc/11.3.0
        export PATH="/scratch.global/saonli/conda-envs/prosper/bin:$PATH"
        export CC="gcc"
        export CXX="g++"
        export CXX11="g++"
        export CXX14="g++"
        export CXX17="g++"
        export CXX20="g++"
        export MAKEFLAGS="CC=gcc CXX=g++ CXX11=g++ CXX14=g++ CXX17=g++ CXX20=g++"
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
          --bfile_tuning {params.tuning_prefix} \
          --pheno_tuning {input.tuning_pheno} \
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
    input:
        resources=rules.preparePRSMethodResources.output.ready,
        env=rules.preparePRSInputs.output.env,
        target_sumstats=rules.preparePROSPERSumstats.output.target_ss,
        training_sumstats=rules.preparePROSPERSumstats.output.training_ss,
        lassosum_target=rules.preparePROSPERLassosumParams.output.target_param,
        lassosum_aux=rules.preparePROSPERLassosumParams.output.aux_param,
        tuning_bed=rules.preparePROSPERTuningTestingSets.output.tuning_bed,
        tuning_bim=rules.preparePROSPERTuningTestingSets.output.tuning_bim,
        tuning_fam=rules.preparePROSPERTuningTestingSets.output.tuning_fam,
        tuning_pheno=rules.preparePROSPERTuningTestingSets.output.tuning_pheno,
        testing_bed=rules.preparePROSPERTuningTestingSets.output.testing_bed,
        testing_bim=rules.preparePROSPERTuningTestingSets.output.testing_bim,
        testing_fam=rules.preparePROSPERTuningTestingSets.output.testing_fam,
        testing_pheno=rules.preparePROSPERTuningTestingSets.output.testing_pheno,
    output:
        done=PRS_METHOD_RUN_DIR / "multi_prosper.done",
    params:
        method="multi_prosper",
        command=prs_method_command_quoted("multi_prosper"),
        extra=prs_method_extra_args("multi_prosper"),
        out_dir=PRS_METHOD_RUN_DIR / "multi_prosper",
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
