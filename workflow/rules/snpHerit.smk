import csv
import json
import re
import subprocess
from pathlib import Path

SNP_HERIT_CONFIG = config.get("snpHerit", {})


def ensure_list(value):
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return list(value)
    return [value]


def run_by_ancestry():
    return bool(SNP_HERIT_CONFIG.get("run_by_ancestry", False))


def active_adjhe_config():
    return bool(
        SNP_HERIT_CONFIG.get("pheno")
        and SNP_HERIT_CONFIG.get("covar")
        and run_by_ancestry()
    )


SNP_HERIT_ACTIVE = active_adjhe_config()

if SNP_HERIT_CONFIG:
    if SNP_HERIT_CONFIG.get("pheno") and not SNP_HERIT_CONFIG.get("covar"):
        raise ValueError(
            "snpHerit.covar must be specified in config when pheno is specified"
        )
    if SNP_HERIT_CONFIG.get("covar") and not SNP_HERIT_CONFIG.get("pheno"):
        raise ValueError(
            "snpHerit.pheno must be specified in config when covar is specified"
        )


if SNP_HERIT_ACTIVE:

    def _out_name():
        out = SNP_HERIT_CONFIG.get("out", "heritability_adjhe.csv")
        return Path(out).name


    def _source_grm_prefix():
        return str(
            SNP_HERIT_CONFIG.get(
                "source_grm_prefix", str(OUT_DIR / "full" / "unrelated")
            )
        )


    def _source_bfile_prefix():
        return str(
            SNP_HERIT_CONFIG.get("source_bfile_prefix", _source_grm_prefix())
        )


    def _plink2_bin():
        return str(
            SNP_HERIT_CONFIG.get(
                "plink2_bin",
                config.get("prsPipeline", {}).get("path_plink2", "plink2"),
            )
        )


    def _ancestry_ref_prefix():
        return str(
            config.get("ancestry", {}).get(
                "reference_panel_prefix",
                REF / "1000G_GRCh38" / "1000G.ensembl.105.with.rsid.gender",
            )
        )


    def _ancestry_ref_pop_file():
        return str(
            config.get("ancestry", {}).get(
                "reference_population_file",
                REF / "1000G_GRCh38" / "1000G.GRCh38.popu",
            )
        )


    def _ancestry_overlap_snplist():
        return str(
            config.get("ancestry", {}).get(
                "overlap_snp_list",
                OUT_DIR
                / "01-globalAncestry"
                / "intermediates"
                / "abcd_1000ggrch38_overlap.snplist",
            )
        )


    def _ancestry_overlap_from_bim():
        return str(config.get("ancestry", {}).get("overlap_from_bim", ""))


    def _ancestry_score_prefix():
        return str(config.get("ancestry", {}).get("score_prefix", "1000ggrch38"))


    def _keep_file_for_subset(wildcards):
        keep = get_ancestry_file(wildcards)
        return keep if keep else []


    # Global ancestry classification should use the high-overlap 1000G GRCh38
    # PLINK panel exposed by the helpers above. The actual PCA/classification
    # rules live outside this file, but these helpers keep the intended
    # reference panel and overlap SNP list explicit so ancestry labels and
    # downstream heritability subsets remain synchronized.


    rule prepareABCDHeritInputsByAncestry:
        container:
            "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
        conda:
            "../../envs/ancNreport.yml"
        threads: 2
        resources:
            nodes=1,
            mem_mb=16000,
            runtime=240,
        input:
            keep=_keep_file_for_subset,
        output:
            pheno_tsv=OUT_DIR
            / "{subset}"
            / "03-snpHeritability"
            / "adjhe_input.pheno.tsv",
            covar_tsv=OUT_DIR
            / "{subset}"
            / "03-snpHeritability"
            / "adjhe_input.covar.tsv",
            pheno_phen=OUT_DIR
            / "{subset}"
            / "03-snpHeritability"
            / "adjhe_input.pheno.phen",
            covar_txt=OUT_DIR
            / "{subset}"
            / "03-snpHeritability"
            / "adjhe_input.covar.txt",
            keep=OUT_DIR / "{subset}" / "03-snpHeritability" / "adjhe_input.keep",
            manifest=OUT_DIR
            / "{subset}"
            / "03-snpHeritability"
            / "adjhe_input_manifest.json",
        params:
            pheno_files=json.dumps(ensure_list(SNP_HERIT_CONFIG.get("pheno"))),
            covar_files=json.dumps(ensure_list(SNP_HERIT_CONFIG.get("covar"))),
            mpheno_names=json.dumps(ensure_list(SNP_HERIT_CONFIG.get("mpheno"))),
            qcovar=json.dumps(ensure_list(SNP_HERIT_CONFIG.get("qcovar"))),
            covar_discrete=json.dumps(
                ensure_list(SNP_HERIT_CONFIG.get("covar_discrete"))
            ),
            iid_col=SNP_HERIT_CONFIG.get("iid_col", "IID"),
            fid_col=SNP_HERIT_CONFIG.get("fid_col", "FID"),
            pheno_filter=json.dumps(SNP_HERIT_CONFIG.get("pheno_filter")),
            covar_filter=json.dumps(SNP_HERIT_CONFIG.get("covar_filter")),
            force_fid_from_iid=bool(SNP_HERIT_CONFIG.get("force_fid_from_iid", True)),
            drop_missing_rows=bool(SNP_HERIT_CONFIG.get("drop_missing_rows", True)),
            source_grm_prefix=_source_grm_prefix(),
        run:
            import pandas as pd

            out_dir = Path(output.pheno_tsv).parent
            out_dir.mkdir(parents=True, exist_ok=True)

            pheno_files = json.loads(params.pheno_files)
            covar_files = json.loads(params.covar_files)
            mpheno_names = json.loads(params.mpheno_names)
            qcovar = json.loads(params.qcovar)
            covar_discrete = json.loads(params.covar_discrete)
            iid_col = params.iid_col
            fid_col = params.fid_col
            pheno_filter = json.loads(params.pheno_filter)
            covar_filter = json.loads(params.covar_filter)
            force_fid_from_iid = params.force_fid_from_iid
            drop_missing_rows = params.drop_missing_rows

            keep_path = str(input.keep) if input.keep else None

            def read_table(path):
                lower = str(path).lower()
                if lower.endswith(".tsv") or lower.endswith(".tab"):
                    return pd.read_csv(path, sep="\t")
                if lower.endswith(".csv"):
                    return pd.read_csv(path)
                if lower.endswith(".parquet"):
                    return pd.read_parquet(path)
                return pd.read_csv(path, sep=r"\s+")

            def apply_optional_filter(df, query_string):
                if not query_string:
                    return df
                try:
                    return df.query(query_string, engine="python")
                except Exception:
                    return df

            def normalize_ids(df):
                cols = {c.lower(): c for c in df.columns}
                iid_key = iid_col.lower()
                fid_key = fid_col.lower()

                if iid_key not in cols:
                    raise ValueError(
                        f"Could not find IID column '{iid_col}' in columns: {list(df.columns)}"
                    )

                df = df.rename(columns={cols[iid_key]: "IID"})

                if fid_key in cols:
                    df = df.rename(columns={cols[fid_key]: "FID"})
                elif "FID" not in df.columns:
                    df["FID"] = df["IID"]

                if force_fid_from_iid:
                    df["FID"] = df["IID"]

                df["FID"] = df["FID"].astype(str)
                df["IID"] = df["IID"].astype(str)
                return df

            def merge_tables(paths, query_string):
                dfs = []
                for path in paths:
                    df = normalize_ids(read_table(path))
                    df = apply_optional_filter(df, query_string)
                    dfs.append(df)

                if not dfs:
                    raise ValueError("No phenotype/covariate files were provided.")

                merged = dfs[0]
                for df in dfs[1:]:
                    keep_cols = [c for c in df.columns if c != "FID"]
                    merged = merged.merge(df[keep_cols], on="IID", how="inner")
                    if "FID_x" in merged.columns:
                        merged["FID"] = merged["FID_x"]
                        drop_cols = [
                            c for c in ["FID_x", "FID_y"] if c in merged.columns
                        ]
                        merged = merged.drop(columns=drop_cols)

                merged["FID"] = merged["IID"]
                merged = merged.drop_duplicates(subset=["IID"])
                return merged

            pheno_df = merge_tables(pheno_files, pheno_filter)
            covar_df = merge_tables(covar_files, covar_filter)

            missing_pheno = [c for c in mpheno_names if c not in pheno_df.columns]
            if missing_pheno:
                raise ValueError(f"Missing phenotype columns: {missing_pheno}")

            all_covar_cols = [c for c in (qcovar + covar_discrete) if c]
            missing_covar = [c for c in all_covar_cols if c not in covar_df.columns]
            if missing_covar:
                raise ValueError(f"Missing covariate columns: {missing_covar}")

            merged = pheno_df[["FID", "IID"] + mpheno_names].merge(
                covar_df[["IID"] + all_covar_cols], on="IID", how="inner"
            )
            merged["FID"] = merged["IID"]

            if keep_path:
                keep_df = pd.read_csv(
                    keep_path, sep=r"\s+", header=None, names=["FID", "IID"], dtype=str
                )
                keep_df["FID"] = keep_df["IID"]
                keep_df = keep_df.drop_duplicates(subset=["FID", "IID"])
                merged = merged.merge(keep_df, on=["FID", "IID"], how="inner")

            grm_id_path = Path(f"{params.source_grm_prefix}.grm.id")
            if not grm_id_path.exists():
                raise ValueError(f"Could not find GRM ID file: {grm_id_path}")

            grm_ids = pd.read_csv(
                grm_id_path,
                sep=r"\s+",
                header=None,
                names=["FID", "IID"],
                dtype=str,
            )
            grm_ids["FID"] = grm_ids["IID"]
            grm_ids["row_order"] = range(len(grm_ids))

            merged = grm_ids[["FID", "IID", "row_order"]].merge(
                merged, on=["FID", "IID"], how="inner"
            )

            before_drop = len(merged)
            if drop_missing_rows:
                merged = merged.dropna(subset=mpheno_names + all_covar_cols)

            merged = merged.sort_values("row_order").reset_index(drop=True)
            merged = merged.drop(columns=["row_order"])

            if merged.empty:
                raise ValueError(
                    f"No samples remain for subset {wildcards.subset} after GRM/phenotype/covariate intersection."
                )

            pheno_out = merged[["FID", "IID"] + mpheno_names].copy()
            covar_out = merged[["FID", "IID"] + all_covar_cols].copy()

            pheno_out.to_csv(output.pheno_tsv, sep="\t", index=False)
            covar_out.to_csv(output.covar_tsv, sep="\t", index=False)
            pheno_out.to_csv(output.pheno_phen, sep=" ", index=False, header=False)
            covar_out.to_csv(output.covar_txt, sep=" ", index=False, header=False)
            merged[["FID", "IID"]].to_csv(
                output.keep, sep=" ", index=False, header=False
            )

            with open(output.manifest, "w") as fh:
                json.dump(
                    {
                        "subset": wildcards.subset,
                        "keep_file": keep_path,
                        "source_grm_prefix": params.source_grm_prefix,
                        "mpheno_names": mpheno_names,
                        "qcovar": qcovar,
                        "covar_discrete": covar_discrete,
                        "n_final": int(len(merged)),
                        "n_before_dropna": int(before_drop),
                    },
                    fh,
                    indent=2,
                )


    rule prepareAdjHEAnalysisSetByAncestry:
        threads: 8
        resources:
            nodes=1,
            mem_mb=64000,
            runtime=1440,
        input:
            keep=rules.prepareABCDHeritInputsByAncestry.output.keep,
            bed=lambda wildcards: f"{_source_bfile_prefix()}.bed",
            bim=lambda wildcards: f"{_source_bfile_prefix()}.bim",
            fam=lambda wildcards: f"{_source_bfile_prefix()}.fam",
        output:
            bed=OUT_DIR
            / "{subset}"
            / "03-snpHeritability"
            / "adjhe_analysis"
            / "unrelated_shared.bed",
            bim=OUT_DIR
            / "{subset}"
            / "03-snpHeritability"
            / "adjhe_analysis"
            / "unrelated_shared.bim",
            fam=OUT_DIR
            / "{subset}"
            / "03-snpHeritability"
            / "adjhe_analysis"
            / "unrelated_shared.fam",
            grm=OUT_DIR
            / "{subset}"
            / "03-snpHeritability"
            / "adjhe_analysis"
            / "unrelated_shared.grm.bin",
            grmid=OUT_DIR
            / "{subset}"
            / "03-snpHeritability"
            / "adjhe_analysis"
            / "unrelated_shared.grm.id",
            grmN=OUT_DIR
            / "{subset}"
            / "03-snpHeritability"
            / "adjhe_analysis"
            / "unrelated_shared.grm.N.bin",
            pca_raw=OUT_DIR
            / "{subset}"
            / "03-snpHeritability"
            / "adjhe_analysis"
            / "unrelated_shared.eigenvec",
            pca_adjhe=OUT_DIR
            / "{subset}"
            / "03-snpHeritability"
            / "adjhe_analysis"
            / "unrelated_shared_adjhe.eigenvec",
        params:
            plink2_bin=_plink2_bin(),
            npc=SNP_HERIT_CONFIG.get("npc", 10),
            source_bfile_prefix=_source_bfile_prefix(),
            prefix=lambda wildcards, output: str(output.bed)[:-4],
            pca_end_col=int(SNP_HERIT_CONFIG.get("npc", 10)) + 2,
        shell:
            """
            set -euo pipefail

            mkdir -p "$(dirname {output.bed})"

            {params.plink2_bin} \
              --bfile {params.source_bfile_prefix} \
              --keep {input.keep} \
              --make-bed \
              --out {params.prefix}

            {params.plink2_bin} \
              --bfile {params.prefix} \
              --make-grm-bin \
              --pca approx {params.npc} \
              --out {params.prefix}

            awk 'NR>1 {{for (i=1; i<={params.pca_end_col}; ++i) printf "%s%s", $i, (i=={params.pca_end_col} ? ORS : OFS)}}' \
              {output.pca_raw} > {output.pca_adjhe}
            """


    rule estimateSnpHeritabilityAdjHEByAncestry:
        threads: 1
        resources:
            nodes=1,
            mem_mb=32000,
            runtime=2880,
        input:
            pheno=rules.prepareABCDHeritInputsByAncestry.output.pheno_phen,
            covar=rules.prepareABCDHeritInputsByAncestry.output.covar_txt,
            manifest=rules.prepareABCDHeritInputsByAncestry.output.manifest,
            grm=rules.prepareAdjHEAnalysisSetByAncestry.output.grm,
            grmid=rules.prepareAdjHEAnalysisSetByAncestry.output.grmid,
            grmN=rules.prepareAdjHEAnalysisSetByAncestry.output.grmN,
            pca=rules.prepareAdjHEAnalysisSetByAncestry.output.pca_adjhe,
        output:
            estimates=OUT_DIR / "{subset}" / "03-snpHeritability" / _out_name(),
        params:
            python_bin=SNP_HERIT_CONFIG.get("python_bin", "python"),
            adjhe_script=SNP_HERIT_CONFIG.get("adjhe_script", ""),
            npc=SNP_HERIT_CONFIG.get("npc", 10),
            grm_prefix=lambda wildcards, input: str(input.grm)[:-8],
        run:
            out_dir = Path(output.estimates).parent
            run_dir = out_dir / "adjhe_runs"
            run_dir.mkdir(parents=True, exist_ok=True)

            if not params.adjhe_script:
                raise ValueError(
                    "snpHerit.adjhe_script must be set for direct AdjHE mode"
                )

            with open(input.manifest) as fh:
                manifest = json.load(fh)

            mpheno_names = manifest["mpheno_names"]
            summary_rows = []

            for idx, pheno_name in enumerate(mpheno_names, start=1):
                safe_name = re.sub(r"[^A-Za-z0-9._-]+", "_", pheno_name)
                per_pheno_out = run_dir / f"{safe_name}.csv"

                cmd = [
                    str(params.python_bin),
                    str(params.adjhe_script),
                    "--prefix",
                    str(params.grm_prefix),
                    "--PC",
                    str(input.pca),
                    "--pheno",
                    str(input.pheno),
                    "--covar",
                    str(input.covar),
                    "--out",
                    str(per_pheno_out),
                    "--npc",
                    str(params.npc),
                    "--mpheno",
                    str(idx),
                ]

                subprocess.run(cmd, check=True)

                with open(per_pheno_out, newline="") as fh:
                    reader = csv.DictReader(fh)
                    rows = list(reader)
                    if not rows:
                        raise ValueError(
                            f"No rows found in AdjHE output: {per_pheno_out}"
                        )
                    for row in rows:
                        row["subset"] = wildcards.subset
                        row["phenotype"] = pheno_name
                        summary_rows.append(row)

            fieldnames = ["subset", "phenotype"]
            for row in summary_rows:
                for key in row.keys():
                    if key not in fieldnames:
                        fieldnames.append(key)

            with open(output.estimates, "w", newline="") as fh:
                writer = csv.DictWriter(fh, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(summary_rows)


    def ancestry_targets():
        ancestries = SNP_HERIT_CONFIG.get("ancestries", ["AFR", "EUR"])
        return expand(
            OUT_DIR / "{anc}" / "03-snpHeritability" / _out_name(),
            anc=ancestries,
    )

    if SNP_HERIT_ACTIVE and run_by_ancestry():

        rule runAncestrySpecificSnpHeritability:
            input:
                rules.classifySamplesByAncestry.output.classifications,
                ancestry_targets()
