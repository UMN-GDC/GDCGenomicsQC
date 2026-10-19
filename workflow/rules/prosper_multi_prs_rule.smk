PROSPER_INPUT_DIR = Path(config.get("prsPipeline", {}).get("generated_input_dir", ""))
PROSPER_RUNTIME_DIR = PROSPER_INPUT_DIR / "prosper" / "PROSPER_runtime"






rule preparePROSPERSumstats:
    log:
        OUT_DIR / "logs" / "preparePROSPERSumstats.log",
    input:
        target_ss=rules.preparePRSInputs.output.target_sumstats,
        training_ss=rules.preparePRSInputs.output.training_sumstats,
        target_bim=rules.preparePRSInputs.output.study_bim,
        training_bim=rules.preparePRSInputs.output.study_anc2_bim,
        reference_bim=Path(
            config.get("prsMethods", {}).get("multi_prosper", {}).get(
                "ld_matrix_dir",
                "/scratch.global/saonli/GDCGenomicsQC/prs_resources/ld/prosper",
            )
        ) / "ref_bim.txt",
    output:
        target_ss=PROSPER_INPUT_DIR / "gwas" / "target_sumstats_PROSPER.txt",
        training_ss=PROSPER_INPUT_DIR / "gwas" / "training_sumstats_PROSPER.txt",
        target_id_map=PROSPER_INPUT_DIR / "metadata" / "AFR_PROSPER_rsid_map.tsv",
    run:
        import pandas as pd

        reference = pd.read_csv(
            str(input.reference_bim),
            sep=r"\s+",
            header=None,
            names=["chr_ref", "rsid_ref", "cm_ref", "pos_ref", "a1_ref", "a0_ref"],
            dtype={"chr_ref": str, "rsid_ref": str, "pos_ref": str},
        )
        reference = reference.drop_duplicates(["chr_ref", "pos_ref"], keep=False)

        complement = str.maketrans("ACGT", "TGCA")

        def make_prosper(ss_file, bim_file, out_file, map_file=None):
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
                names=["chr_bim", "rsid_orig", "cm", "pos_bim", "a1_bim", "a0_bim"],
                dtype={"chr_bim": str, "rsid_orig": str, "pos_bim": str},
            )
            bim = bim.drop_duplicates("rsid_orig", keep=False)

            required = ["rsid", "a1", "a0", "beta", "beta_se", "n_eff"]
            missing = [c for c in required if c not in ss.columns]
            if missing:
                raise ValueError(f"Missing columns in {ss_file}: {', '.join(missing)}")

            ss = ss.merge(
                bim[["chr_bim", "rsid_orig", "pos_bim"]],
                left_on="rsid",
                right_on="rsid_orig",
                how="inner",
            )
            ss = ss.merge(
                reference,
                left_on=["chr_bim", "pos_bim"],
                right_on=["chr_ref", "pos_ref"],
                how="inner",
            )
            for allele in ("a1", "a0", "a1_ref", "a0_ref"):
                ss[allele] = ss[allele].astype(str).str.upper()

            same = (ss["a1"] == ss["a1_ref"]) & (ss["a0"] == ss["a0_ref"])
            swapped = (ss["a1"] == ss["a0_ref"]) & (ss["a0"] == ss["a1_ref"])
            comp_a1 = ss["a1"].str.translate(complement)
            comp_a0 = ss["a0"].str.translate(complement)
            comp_same = (comp_a1 == ss["a1_ref"]) & (comp_a0 == ss["a0_ref"])
            comp_swapped = (comp_a1 == ss["a0_ref"]) & (comp_a0 == ss["a1_ref"])
            keep = same | swapped | comp_same | comp_swapped
            ss = ss.loc[keep].copy()
            flipped = swapped | comp_swapped
            ss.loc[flipped.loc[ss.index], "beta"] *= -1

            if map_file is not None:
                rsid_map = ss[["rsid_orig", "rsid_ref"]].drop_duplicates("rsid_orig", keep=False)
                rsid_map.to_csv(str(map_file), sep="\t", index=False, header=False)

            result = ss.loc[
                ss["beta"].notna() & ss["beta_se"].notna(),
                ["rsid_ref", "chr_ref", "a1_ref", "a0_ref", "beta", "beta_se", "n_eff"],
            ].copy()
            result.columns = ["rsid", "chr", "a1", "a0", "beta", "beta_se", "n_eff"]
            result.to_csv(str(out_file), sep="\t", index=False)
            return len(result)

        target_n = make_prosper(
            input.target_ss, input.target_bim, output.target_ss, output.target_id_map
        )
        training_n = make_prosper(input.training_ss, input.training_bim, output.training_ss)

        with open(str(log), "w", encoding="utf-8") as handle:
            handle.write(f"Mapped {target_n} target SNPs to PROSPER rsIDs.\n")
            handle.write(f"Mapped {training_n} training SNPs to PROSPER rsIDs.\n")



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
        shared_tuning_keep=rules.preparePRSInputs.output.target_tuning_keep,
        shared_testing_keep=rules.preparePRSInputs.output.target_testing_keep,
        shared_tuning_pheno=rules.preparePRSInputs.output.target_tuning_pheno,
        shared_testing_pheno=rules.preparePRSInputs.output.target_testing_pheno,
        target_id_map=rules.preparePROSPERSumstats.output.target_id_map,
        reference_bim=Path(
            config.get("prsMethods", {}).get("multi_prosper", {}).get(
                "ld_matrix_dir",
                "/scratch.global/saonli/GDCGenomicsQC/prs_resources/ld/prosper",
            )
        ) / "ref_bim.txt",
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
        runtime_package=directory(PROSPER_RUNTIME_DIR),
    params:
        plink2=config.get("prsMethods", {}).get("plink2", "/projects/standard/gdc/public/plink2"),
        source_package=config.get("prsMethods", {}).get("multi_prosper", {}).get(
            "software_dir",
            "/scratch.global/saonli/GDCGenomicsQC/prs_resources/software/PROSPER",
        ),
        ld_ref_dir=config.get("prsMethods", {}).get("multi_prosper", {}).get(
            "ld_matrix_dir",
            "/scratch.global/saonli/GDCGenomicsQC/prs_resources/ld/prosper",
        ),
        study_prefix=lambda wc, input: str(input.bed)[:-4],
        tuning_prefix=lambda wc, output: str(output.tuning_bed)[:-4],
        testing_prefix=lambda wc, output: str(output.testing_bed)[:-4],
    shell:
        """
        set -euo pipefail

        mkdir -p $(dirname {output.tuning_keep})
        mkdir -p $(dirname {output.tuning_bed})
        rm -rf {output.runtime_package}
        cp -a {params.source_package} {output.runtime_package}
        rm -rf {output.runtime_package}/AFR {output.runtime_package}/EUR
        ln -s {params.ld_ref_dir}/AFR {output.runtime_package}/AFR
        ln -s {params.ld_ref_dir}/EUR {output.runtime_package}/EUR
        cp {input.reference_bim} {output.runtime_package}/ref_bim.txt

        # PROSPER's lassosum2 script uses `%dopar%` even when NCORES=1.  On
        # MSI this hides chromosome errors inside a foreach worker, so run the
        # disposable analysis copy serially instead of modifying shared code.
        grep -q '%dopar%' {output.runtime_package}/scripts/lassosum2.R
        sed -i 's/%dopar%/%do%/g' {output.runtime_package}/scripts/lassosum2.R

        # `%do%` still wraps errors as "task 1 failed".  Make the fitting loop
        # a regular R for-loop, so any invalid LD block is reported precisely.
        grep -q 'ff <- foreach(j = 1:length(allchrom)' {output.runtime_package}/scripts/lassosum2.R
        perl -0pi -e 's/\Q  ff <- foreach(j = 1:length(allchrom), ii = icount(), .final = function(x) NULL) %do% {{\E\n/  for (j in seq_along(allchrom)) {{\n    ii <- j\n/' {output.runtime_package}/scripts/lassosum2.R
        grep -q 'for (j in seq_along(allchrom))' {output.runtime_package}/scripts/lassosum2.R

        # PLINK 2 writes five metadata columns before multi-score results.
        # The upstream script drops four, incorrectly treating the dosage total
        # as a 26th PRS and selecting an out-of-range weight column.
        grep -q 'SCORE <- SCORE\[,-1:-4\]' {output.runtime_package}/scripts/lassosum2.R
        sed -i 's/SCORE <- SCORE\[,-1:-4\]/SCORE <- SCORE[,-1:-5]/' {output.runtime_package}/scripts/lassosum2.R

        # Reuse the outcome-blinded split shared by every benchmark method.
        cp {input.shared_tuning_keep} {output.tuning_keep}
        cp {input.shared_testing_keep} {output.testing_keep}
        # PROSPER's upstream scripts expect a headerless PLINK phenotype file.
        awk 'NR > 1 {{print $1, $2, $3}}' {input.shared_tuning_pheno} > {output.tuning_pheno}
        awk 'NR > 1 {{print $1, $2, $3}}' {input.shared_testing_pheno} > {output.testing_pheno}

        {params.plink2} \
          --bfile {params.study_prefix} \
          --keep {output.tuning_keep} \
          --update-name {input.target_id_map} 2 1 \
          --make-bed \
          --out {params.tuning_prefix}

        {params.plink2} \
          --bfile {params.study_prefix} \
          --keep {output.testing_keep} \
          --update-name {input.target_id_map} 2 1 \
          --make-bed \
          --out {params.testing_prefix}

        echo "Prepared PROSPER tuning/testing PLINK sets" > {log}
        echo "Using LD panels from: {params.ld_ref_dir}" >> {log}
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
        runtime_package=rules.preparePROSPERTuningTestingSets.output.runtime_package,
    output:
        target_param=PRS_METHOD_RUN_DIR / "multi_prosper" / "lassosum2_AFR" / "param.txt",
        aux_param=PRS_METHOD_RUN_DIR / "multi_prosper" / "lassosum2_EUR" / "param.txt",
    params:
        package=lambda wc, input: str(input.runtime_package),
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
        runtime_package=rules.preparePROSPERTuningTestingSets.output.runtime_package,
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
