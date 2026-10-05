#!/usr/bin/env python3
"""Wrapper around PredLMM-ACE: align GRM/knots/metadata and run predlmm-fit."""
import argparse, json, re, subprocess, sys
from pathlib import Path
import pandas as pd


def read_table(path):
    p = Path(path)
    ext = p.suffix.lower()
    if ext in (".csv",):
        return pd.read_csv(p, dtype=str)
    if ext in (".tsv", ".tab"):
        return pd.read_csv(p, sep="\t", dtype=str)
    return pd.read_csv(p, sep=r"\s+", dtype=str)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grm-prefix", required=True)
    ap.add_argument("--pheno", required=True)
    ap.add_argument("--pheno-col", default=None, help="phenotype column name; defaults to first column after FID/IID")
    ap.add_argument("--covar", default=None)
    ap.add_argument("--eigenvec", default=None)
    ap.add_argument("--npc", type=int, default=10)
    ap.add_argument("--family-column", default="FID")
    ap.add_argument("--model", default="auto", choices=["auto", "ACE", "AE"])
    ap.add_argument("--rank", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=20260918)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    out = Path(args.out)
    work = out.parent / (out.stem + "_predlm_ace")
    work.mkdir(parents=True, exist_ok=True)

    pheno = read_table(args.pheno)
    if args.pheno_col is None:
        args.pheno_col = pheno.columns[2]
    need = {"FID", "IID", args.pheno_col}
    missing = need - set(pheno.columns)
    if missing:
        sys.exit(f"phenotype file missing columns: {missing}")
    bits = [pheno[["FID", "IID", args.pheno_col]]]

    cont_extra, cat_extra = [], []
    if args.covar:
        covar = read_table(args.covar)
        covar = covar.rename(columns={c: c for c in covar.columns})
        for c in covar.columns:
            if c in ("FID", "IID"):
                continue
            try:
                pd.to_numeric(covar[c])
                cont_extra.append(c)
            except (ValueError, TypeError):
                cat_extra.append(c)
        bits.append(covar)

    eigen_cols = []
    if args.eigenvec:
        ev = read_table(args.eigenvec)
        first = ev.columns[0]
        if first.startswith("#"):
            ev = ev.rename(columns={first: first.lstrip("#")})
        pc_cols = [c for c in ev.columns if c.upper().startswith("PC")][: args.npc]
        eigen_cols = pc_cols
        bits.append(ev[["FID", "IID"] + pc_cols])

    if args.family_column != "FID" and args.family_column not in set().union(*[set(b.columns) for b in bits]):
        sys.exit(f"family column {args.family_column} not found in phenotype/covariate files")

    meta = bits[0]
    for b in bits[1:]:
        meta = meta.merge(b, on=["FID", "IID"], how="inner")
    used = [args.pheno_col, args.family_column if args.family_column in meta.columns else "FID"] + cont_extra + cat_extra + eigen_cols
    meta = meta.dropna(subset=used)
    if meta.empty:
        sys.exit("no samples with complete phenotype/covariate/PC data")

    # ensure family column exists (family_column == 'FID' is already present)
    meta_path = work / "metadata.tsv"
    meta.to_csv(meta_path, sep="\t", index=False)

    knots = work / "knots.tsv"
    subprocess.run(["predlmm-select-knots",
        "--metadata", str(meta_path), "--grm-id", args.grm_prefix + ".grm.id",
        "--family-column", args.family_column, "--rank", str(args.rank),
        "--seed", str(args.seed), "--out", str(knots)], check=True)

    factor_dir = work / "factor"
    subprocess.run(["predlmm-grm-to-nystrom", "--grm-prefix", args.grm_prefix,
                    "--knots", str(knots), "--out", str(factor_dir)], check=True)

    continuous = (eigen_cols + cont_extra)
    fit_out = work / "fit.json"
    def run_fit(model):
        cmd = ["predlmm-fit", "--factor-dir", str(factor_dir), "--metadata", str(meta_path),
               "--phenotype", args.pheno_col, "--family-column", args.family_column,
               "--model", model, "--standard-errors", "--out", str(fit_out)]
        if continuous:
            cmd += ["--continuous", ",".join(continuous)]
        if cat_extra:
            cmd += ["--categorical", ",".join(cat_extra)]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        return proc

    chosen = args.model
    proc = run_fit("ACE") if chosen in ("auto", "ACE") else run_fit("AE")
    if chosen == "auto" and proc.returncode != 0 and re.search(r"not identifiable|singleton", proc.stderr + proc.stdout):
        chosen = "AE"
        proc = run_fit("AE")
    if proc.returncode != 0:
        sys.exit(f"predlmm-fit failed:\n{proc.stdout}\n{proc.stderr}")

    res = json.loads(fit_out.read_text())
    row = {k: res.get(k) for k in ["model", "n", "rank", "h2", "c2", "se_h2", "se_c2",
                                   "converged", "iterations", "message", "phenotype",
                                   "family_column", "n_families", "n_multi_participant_families"]}
    row["model_requested"] = args.model
    row["model_used"] = chosen
    df = pd.DataFrame([row])
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
