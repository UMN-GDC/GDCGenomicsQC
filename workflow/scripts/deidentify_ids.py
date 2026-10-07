#!/usr/bin/env python3
"""
De-identify files using a crosswalk (raw_id -> release_id + suffix).

Input file types supported:
- Text/CSV/TSV with a sample_id column containing raw IDs (e.g., "pscid+C" or "array_channel_pscidC")
- PLINK .fam/.psam (FID/IID columns)
- GRM .grm.id (FID IID)
- Eigenvec (FID IID PC1...)
- CNV files (sample_id column with format like "array_channel_pscidC")

The crosswalk maps pscid -> (release_candid, suffix).
Sample IDs in input files may have various formats; this extracts pscid + suffix, maps pscid,
and rebuilds the ID as release_candid + suffix.
"""
import argparse, sys, re, pandas as pd
from pathlib import Path


def read_table(path, **kwargs):
    p = Path(path)
    ext = p.suffix.lower()
    if ext in (".csv",):
        return pd.read_csv(p, **kwargs)
    if ext in (".tsv", ".tab"):
        return pd.read_csv(p, sep="\t", **kwargs)
    if ext in (".xlsx", ".xls"):
        return pd.read_excel(p, **kwargs)
    return pd.read_csv(p, sep=r"\s+", **kwargs)


def write_table(df, path, orig_path):
    p = Path(path)
    ext = p.suffix.lower()
    if ext in (".csv",):
        df.to_csv(p, index=False)
    elif ext in (".tsv", ".tab"):
        df.to_csv(p, sep="\t", index=False)
    else:
        # Default: space-separated, no header for PLINK-style
        df.to_csv(p, sep=" ", index=False, header=False)


def extract_pscid_suffix(sample_id, pattern=None):
    """
    Extract (pscid, suffix) from a sample ID string.
    Default pattern: digits followed by C or M at end.
    If pattern provided, use that regex with groups (pscid, suffix).
    """
    s = str(sample_id).strip()
    if pattern:
        m = re.match(pattern, s)
        if m:
            return m.group(1), m.group(2)
    # Default: look for 8+ digits followed by C or M at end
    m = re.search(r"(\d{8,})([CM])$", s)
    if m:
        return m.group(1), m.group(2)
    # Fallback: split on last underscore or just return as-is
    return s, "C"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="Input file to de-identify")
    ap.add_argument("--output", required=True, help="Output file")
    ap.add_argument("--crosswalk", required=True, help="Crosswalk: pscid -> release_candid (+ optional suffix)")
    ap.add_argument("--pscid-col", default="pscid", help="Column name for pscid in crosswalk")
    ap.add_argument("--release-col", default="release_candid", help="Column name for release_candid in crosswalk")
    ap.add_argument("--suffix-col", default=None, help="Optional suffix column in crosswalk")
    ap.add_argument("--sample-id-col", default=None, help="Column containing sample IDs in input (default: auto-detect)")
    ap.add_argument("--id-pattern", default=None, help="Regex pattern with 2 groups (pscid, suffix) for sample IDs")
    ap.add_argument("--file-type", choices=["auto", "fam", "psam", "grm_id", "eigenvec", "cnv", "generic"], default="auto",
                    help="File type for column handling")
    args = ap.parse_args()

    # Load crosswalk
    cw = pd.read_csv(args.crosswalk, dtype=str)
    if args.pscid_col not in cw.columns:
        sys.exit(f"pscid column '{args.pscid_col}' not in crosswalk. Columns: {list(cw.columns)}")
    if args.release_col not in cw.columns:
        sys.exit(f"release column '{args.release_col}' not in crosswalk. Columns: {list(cw.columns)}")

    cw[args.pscid_col] = cw[args.pscid_col].astype(str)
    cw[args.release_col] = cw[args.release_col].astype(str)

    if args.suffix_col and args.suffix_col in cw.columns:
        cw[args.suffix_col] = cw[args.suffix_col].astype(str)
        mapping = dict(zip(cw[args.pscid_col], zip(cw[args.release_col], cw[args.suffix_col])))
    else:
        mapping = {row[args.pscid_col]: (row[args.release_col], "C") for _, row in cw.iterrows()}

    print(f"Loaded crosswalk: {len(mapping)} mappings")

    # Load input
    df = read_table(args.input, dtype=str)
    print(f"Loaded input: {df.shape[0]} rows, {df.shape[1]} cols")

    # Determine file type and ID columns
    if args.file_type == "auto":
        cols = set(df.columns.str.lower())
        if {"fid", "iid"}.issubset(cols):
            args.file_type = "fam"
        elif {"iid", "#iid"}.intersection(cols):
            args.file_type = "psam"
        elif set(df.columns[:2].str.lower()) == {"fid", "iid"}:
            args.file_type = "grm_id"
        elif "sample_id" in df.columns:
            args.file_type = "cnv"
        elif "fid" in cols and "iid" in cols:
            args.file_type = "eigenvec"
        else:
            args.file_type = "generic"

    # Process based on file type
    def deid_id(raw):
        pscid, suf = extract_pscid_suffix(raw, args.id_pattern)
        if pscid in mapping:
            rel, mapped_suf = mapping[pscid]
            return rel + (mapped_suf if mapped_suf else suf)
        # Fallback: keep original if not in crosswalk
        return raw

    if args.file_type in ("fam", "psam", "grm_id", "eigenvec"):
        # FID/IID columns
        fid_col = "FID" if "FID" in df.columns else ("#FID" if "#FID" in df.columns else df.columns[0])
        iid_col = "IID" if "IID" in df.columns else ("#IID" if "#IID" in df.columns else df.columns[1])

        df[fid_col] = df[fid_col].apply(deid_id)
        df[iid_col] = df[iid_col].apply(deid_id)
        print(f"De-identified {fid_col}/{iid_col}")

    elif args.file_type == "cnv":
        # sample_id column
        sid_col = args.sample_id_col or "sample_id"
        if sid_col not in df.columns:
            sys.exit(f"sample_id column '{sid_col}' not found. Columns: {list(df.columns)}")
        df[sid_col] = df[sid_col].apply(deid_id)
        print(f"De-identified {sid_col}")

    else:
        # Generic: try to find ID column
        sid_col = args.sample_id_col
        if not sid_col:
            for c in ["sample_id", "IID", "iid", "ID", "id"]:
                if c in df.columns:
                    sid_col = c
                    break
        if sid_col and sid_col in df.columns:
            df[sid_col] = df[sid_col].apply(deid_id)
            print(f"De-identified {sid_col}")
        else:
            sys.exit(f"Could not identify ID column for generic file. Columns: {list(df.columns)}")

    # Write output
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    write_table(df, args.output, args.input)
    print(f"Wrote de-identified file: {args.output}")


if __name__ == "__main__":
    main()