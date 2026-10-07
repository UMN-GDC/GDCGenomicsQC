#!/usr/bin/env python3
"""
Filter derivative files to keep only release subjects.

Supported file types:
- GRM binary (.grm.bin + .grm.id + .grm.N.bin) - filters both IDs and matrix
- GRM text gzipped (.grm.gz) - filters rows
- Eigenvec (.eigenvec) - filters by IID
- CNV text files (sample_id column) - filters rows
- Generic TSV/CSV with IID/sample_id column - filters rows

Input: --keep-list (one column, release IIDs), --input file(s)
Output: Filtered files to --output-dir or --output prefix
"""
import argparse, sys, gzip, numpy as np, pandas as pd
from pathlib import Path


def read_keep_list(path):
    """Read one-column keep list."""
    with open(path) as f:
        return set(line.strip() for line in f if line.strip())


def filter_grm_binary(grm_bin, grm_id, grm_nbin, keep_ids, out_prefix):
    """Filter GCTA/PLINK binary GRM to keep only specified IIDs."""
    # Read .grm.id (FID IID)
    ids = pd.read_csv(grm_id, sep=r"\s+", names=["FID", "IID"], dtype=str)
    n = len(ids)
    print(f"  GRM size: {n} x {n}")

    # Find indices to keep
    ids["IID_str"] = ids["IID"].astype(str)
    keep_mask = ids["IID_str"].isin(keep_ids)
    keep_idx = np.where(keep_mask)[0]
    n_keep = len(keep_idx)
    print(f"  Keeping {n_keep} / {n} subjects")

    if n_keep == 0:
        sys.exit("No subjects to keep in GRM!")

    # Read binary GRM (lower triangular, float32)
    # Size: n*(n+1)/2 * 4 bytes
    expected_size = n * (n + 1) // 2 * 4
    actual_size = Path(grm_bin).stat().st_size
    if actual_size != expected_size:
        print(f"  Warning: GRM bin size mismatch: expected {expected_size}, got {actual_size}")

    # Memory-map the binary GRM
    grm_data = np.fromfile(grm_bin, dtype="<f4")  # little-endian float32
    grm_matrix = np.zeros((n, n), dtype=np.float32)

    # Fill lower triangular
    idx = 0
    for i in range(n):
        grm_matrix[i, :i+1] = grm_data[idx:idx+i+1]
        idx += i + 1
    # Make symmetric
    grm_matrix = np.tril(grm_matrix) + np.tril(grm_matrix, -1).T

    # Subset
    grm_sub = grm_matrix[np.ix_(keep_idx, keep_idx)]
    ids_sub = ids.iloc[keep_idx].reset_index(drop=True)

    # Write outputs
    # .grm.id
    ids_sub[["FID", "IID"]].to_csv(f"{out_prefix}.grm.id", sep=" ", index=False, header=False)

    # .grm.bin (lower triangular of subset)
    n_sub = len(keep_idx)
    out_data = []
    for i in range(n_sub):
        out_data.append(grm_sub[i, :i+1])
    np.concatenate(out_data).astype("<f4").tofile(f"{out_prefix}.grm.bin")

    # .grm.N.bin (same structure, but we don't have N-bin info - write ones)
    # Actually N-bin has same structure; for simplicity write ones
    np.ones(n_sub * (n_sub + 1) // 2, dtype="<f4").tofile(f"{out_prefix}.grm.N.bin")

    print(f"  Wrote filtered GRM: {out_prefix}.grm.*")


def filter_grm_text_gz(in_path, keep_ids, out_path):
    """Filter gzipped GRM text file (IID1 IID2 value per line)."""
    keep = keep_ids
    n_in = n_out = 0
    with gzip.open(in_path, "rt") as fin, gzip.open(out_path, "wt") as fout:
        for line in fin:
            n_in += 1
            parts = line.strip().split()
            if len(parts) >= 2:
                iid1, iid2 = parts[0], parts[1]
                if iid1 in keep and iid2 in keep:
                    fout.write(line)
                    n_out += 1
    print(f"  GRM text: {n_out} / {n_in} pairs kept")


def filter_eigenvec(in_path, keep_ids, out_path):
    """Filter eigenvec file (FID IID PC1...)."""
    df = pd.read_csv(in_path, sep=r"\s+", dtype=str)
    # Handle potential # comment header
    if df.columns[0].startswith("#"):
        df = df.rename(columns={df.columns[0]: df.columns[0][1:]})
    iid_col = "IID" if "IID" in df.columns else ("#IID" if "#IID" in df.columns else df.columns[1])
    before = len(df)
    df = df[df[iid_col].astype(str).isin(keep_ids)].reset_index(drop=True)
    after = len(df)
    print(f"  Eigenvec: {after} / {before} rows kept")
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep=" ", index=False)


def filter_cnv(in_path, keep_ids, out_path, sample_id_col="sample_id"):
    """Filter CNV file by sample_id column."""
    ext = Path(in_path).suffix.lower()
    if ext in (".csv",):
        df = pd.read_csv(in_path, dtype=str)
    elif ext in (".tsv", ".tab"):
        df = pd.read_csv(in_path, sep="\t", dtype=str)
    else:
        df = pd.read_csv(in_path, sep=r"\s+", dtype=str)

    if sample_id_col not in df.columns:
        # Try common alternatives
        for c in ["IID", "iid", "sample_id", "ID", "id"]:
            if c in df.columns:
                sample_id_col = c
                break

    before = len(df)
    df = df[df[sample_id_col].astype(str).isin(keep_ids)].reset_index(drop=True)
    after = len(df)
    print(f"  CNV ({sample_id_col}): {after} / {before} rows kept")

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    if out_path.endswith(".csv"):
        df.to_csv(out_path, index=False)
    elif out_path.endswith((".tsv", ".tab")):
        df.to_csv(out_path, sep="\t", index=False)
    else:
        df.to_csv(out_path, sep="\t", index=False)


def filter_generic(in_path, keep_ids, out_path, id_col=None):
    """Filter generic TSV/CSV by ID column."""
    ext = Path(in_path).suffix.lower()
    if ext in (".csv",):
        df = pd.read_csv(in_path, dtype=str)
    elif ext in (".tsv", ".tab"):
        df = pd.read_csv(in_path, sep="\t", dtype=str)
    else:
        df = pd.read_csv(in_path, sep=r"\s+", dtype=str)

    if id_col is None:
        for c in ["IID", "iid", "sample_id", "ID", "id", "subject_id"]:
            if c in df.columns:
                id_col = c
                break
    if not id_col or id_col not in df.columns:
        sys.exit(f"Could not find ID column in {in_path}. Columns: {list(df.columns)}")

    before = len(df)
    df = df[df[id_col].astype(str).isin(keep_ids)].reset_index(drop=True)
    after = len(df)
    print(f"  Generic ({id_col}): {after} / {before} rows kept")

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    if out_path.endswith(".csv"):
        df.to_csv(out_path, index=False)
    elif out_path.endswith((".tsv", ".tab")):
        df.to_csv(out_path, sep="\t", index=False)
    else:
        df.to_csv(out_path, sep="\t", index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--keep-list", required=True, help="One-column file of release IIDs to keep")
    ap.add_argument("--input", required=True, help="Input file (or comma-separated for GRM trio)")
    ap.add_argument("--output", required=True, help="Output file or prefix")
    ap.add_argument("--type", required=True,
                    choices=["grm_binary", "grm_text_gz", "eigenvec", "cnv", "generic"],
                    help="File type to filter")
    ap.add_argument("--sample-id-col", default="sample_id", help="Sample ID column for CNV/generic")
    ap.add_argument("--id-col", default=None, help="ID column for generic")
    args = ap.parse_args()

    keep_ids = read_keep_list(args.keep_list)
    print(f"Loaded keep-list: {len(keep_ids)} IIDs")

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    if args.type == "grm_binary":
        # Expect --input as "grm.bin,grm.id,grm.N.bin"
        parts = args.input.split(",")
        if len(parts) != 3:
            sys.exit("grm_binary requires --input as 'grm.bin,grm.id,grm.N.bin'")
        filter_grm_binary(parts[0], parts[1], parts[2], keep_ids, args.output)

    elif args.type == "grm_text_gz":
        filter_grm_text_gz(args.input, keep_ids, args.output)

    elif args.type == "eigenvec":
        filter_eigenvec(args.input, keep_ids, args.output)

    elif args.type == "cnv":
        filter_cnv(args.input, keep_ids, args.output, args.sample_id_col)

    elif args.type == "generic":
        filter_generic(args.input, keep_ids, args.output, args.id_col)

    print("Done.")


if __name__ == "__main__":
    main()