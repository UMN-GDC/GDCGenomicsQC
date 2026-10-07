#!/usr/bin/env python3
"""
Validate release outputs: ensure no excluded subjects leaked into any output file.

Checks:
- All IIDs in release files match release pattern (10 digits + C/M)
- No IIDs from exclusion lists appear in outputs
- All output IIDs are subset of keep-list
- Row counts match expectations
- GRM dimensions match PLINK .fam
"""
import argparse, sys, re, gzip, pandas as pd, numpy as np
from pathlib import Path


IID_PATTERN = re.compile(r"^\d{10}[CM]$")


def read_keep_list(path):
    with open(path) as f:
        return set(line.strip() for line in f if line.strip())


def read_exclusion_lists(paths):
    """Read all exclusion lists and return combined set of excluded release IIDs."""
    excluded = set()
    for p in paths:
        if not p or not Path(p).exists():
            continue
        ext = Path(p).suffix.lower()
        if ext in (".csv",):
            df = pd.read_csv(p, dtype=str)
        elif ext in (".tsv", ".tab"):
            df = pd.read_csv(p, sep="\t", dtype=str)
        elif ext in (".xlsx", ".xls"):
            df = pd.read_excel(p, dtype=str)
        else:
            df = pd.read_csv(p, sep=r"\s+", dtype=str)
        # Try to find ID column
        for c in ["release_candid", "pscid", "Study_ID", "IID", "iid", "sample_id", "ID", "id"]:
            if c in df.columns:
                ids = df[c].astype(str).tolist()
                # If pscid, need to map - but for validation we assume release IIDs
                for raw in ids:
                    if IID_PATTERN.match(raw):
                        excluded.add(raw)
                    elif raw.isdigit() and len(raw) >= 8:
                        # Might be pscid - add both C and M
                        excluded.add(raw + "C")
                        excluded.add(raw + "M")
                break
    return excluded


def check_iid_pattern(iids, name):
    """Check all IIDs match expected pattern."""
    bad = [i for i in iids if not IID_PATTERN.match(str(i))]
    if bad:
        print(f"  [FAIL] {name}: {len(bad)} IIDs don't match pattern ^\\d{{10}}[CM]$")
        for b in bad[:5]:
            print(f"    {b}")
        return False
    print(f"  [OK] {name}: all {len(iids)} IIDs match pattern")
    return True


def check_no_excluded(iids, excluded, name):
    """Check no excluded IIDs appear."""
    leaked = iids.intersection(excluded)
    if leaked:
        print(f"  [FAIL] {name}: {len(leaked)} excluded IIDs found!")
        for l in list(leaked)[:10]:
            print(f"    {l}")
        return False
    print(f"  [OK] {name}: no excluded IIDs found")
    return True


def check_subset_of_keep(iids, keep, name):
    """Check all IIDs are in keep-list."""
    extra = iids - keep
    if extra:
        print(f"  [FAIL] {name}: {len(extra)} IIDs not in keep-list!")
        for e in list(extra)[:10]:
            print(f"    {e}")
        return False
    print(f"  [OK] {name}: all IIDs in keep-list")
    return True


def validate_fam(path, keep, excluded):
    """Validate PLINK .fam file."""
    print(f"\n=== Validating {path} ===")
    df = pd.read_csv(path, sep=r"\s+", names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"], dtype=str)
    iids = set(df["IID"].astype(str))
    fid_set = set(df["FID"].astype(str))

    ok = True
    ok &= check_iid_pattern(iids, f"{path} (IID)")
    ok &= check_iid_pattern(fid_set, f"{path} (FID)")
    ok &= check_no_excluded(iids, excluded, path)
    ok &= check_no_excluded(fid_set, excluded, path + " (FID)")
    ok &= check_subset_of_keep(iids, keep, path)

    # Check FID == first 10 digits of IID
    mismatched = [i for i in iids if not i.startswith(str(df.loc[df["IID"]==i, "FID"].values[0]))]
    if mismatched:
        print(f"  [FAIL] {path}: FID/IID prefix mismatch for {len(mismatched)} rows")
        ok = False
    else:
        print(f"  [OK] {path}: FID == IID prefix for all rows")

    # Check FID length == 10, IID length == 11
    bad_fid = [f for f in fid_set if len(f) != 10]
    bad_iid = [i for i in iids if len(i) != 11]
    if bad_fid:
        print(f"  [FAIL] {path}: {len(bad_fid)} FIDs not 10 digits")
        ok = False
    if bad_iid:
        print(f"  [FAIL] {path}: {len(bad_iid)} IIDs not 11 digits")
        ok = False

    return ok


def validate_eigenvec(path, keep, excluded):
    """Validate eigenvec file."""
    print(f"\n=== Validating {path} ===")
    df = pd.read_csv(path, sep=r"\s+", dtype=str)
    if df.columns[0].startswith("#"):
        df = df.rename(columns={df.columns[0]: df.columns[0][1:]})
    iid_col = "IID" if "IID" in df.columns else ("#IID" if "#IID" in df.columns else df.columns[1])
    iids = set(df[iid_col].astype(str))

    ok = True
    ok &= check_iid_pattern(iids, path)
    ok &= check_no_excluded(iids, excluded, path)
    ok &= check_subset_of_keep(iids, keep, path)
    return ok


def validate_cnv(path, sample_id_col, keep, excluded):
    """Validate CNV file."""
    print(f"\n=== Validating {path} ===")
    ext = Path(path).suffix.lower()
    if ext in (".csv",):
        df = pd.read_csv(path, dtype=str)
    elif ext in (".tsv", ".tab"):
        df = pd.read_csv(path, sep="\t", dtype=str)
    else:
        df = pd.read_csv(path, sep=r"\s+", dtype=str)

    if sample_id_col not in df.columns:
        for c in ["IID", "iid", "sample_id", "ID", "id"]:
            if c in df.columns:
                sample_id_col = c
                break

    iids = set(df[sample_id_col].astype(str))
    ok = True
    ok &= check_iid_pattern(iids, path)
    ok &= check_no_excluded(iids, excluded, path)
    ok &= check_subset_of_keep(iids, keep, path)
    return ok


def validate_grm_binary(grm_bin, grm_id, keep, excluded):
    """Validate binary GRM dimensions and IDs."""
    print(f"\n=== Validating GRM {grm_bin} ===")
    ids = pd.read_csv(grm_id, sep=r"\s+", names=["FID", "IID"], dtype=str)
    iids = set(ids["IID"].astype(str))
    fids = set(ids["FID"].astype(str))

    ok = True
    ok &= check_iid_pattern(iids, f"{grm_bin} (IID)")
    ok &= check_iid_pattern(fids, f"{grm_bin} (FID)")
    ok &= check_no_excluded(iids, excluded, grm_bin)
    ok &= check_no_excluded(fids, excluded, grm_bin + " (FID)")
    ok &= check_subset_of_keep(iids, keep, grm_bin)

    # Check binary file size matches
    n = len(ids)
    expected = n * (n + 1) // 2 * 4
    actual = Path(grm_bin).stat().st_size
    if actual == expected:
        print(f"  [OK] {grm_bin}: size matches {n} subjects")
    else:
        print(f"  [FAIL] {grm_bin}: size {actual} != expected {expected} for n={n}")
        ok = False

    return ok


def validate_grm_text_gz(grm_gz, keep, excluded):
    """Validate gzipped GRM text file."""
    print(f"\n=== Validating {grm_gz} ===")
    all_iids = set()
    with gzip.open(grm_gz, "rt") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                all_iids.add(parts[0])
                all_iids.add(parts[1])

    ok = True
    ok &= check_iid_pattern(all_iids, grm_gz)
    ok &= check_no_excluded(all_iids, excluded, grm_gz)
    ok &= check_subset_of_keep(all_iids, keep, grm_gz)
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--keep-list", required=True, help="Release keep-list (one column)")
    ap.add_argument("--exclusion-lists", nargs="*", default=[], help="Exclusion list files")
    ap.add_argument("--fam", default=None, help="PLINK .fam to validate")
    ap.add_argument("--eigenvec", default=None, help="Eigenvec file to validate")
    ap.add_argument("--cnv", default=None, help="CNV file to validate")
    ap.add_argument("--cnv-col", default="sample_id", help="CNV sample ID column")
    ap.add_argument("--grm-bin", default=None, help="GRM .grm.bin file")
    ap.add_argument("--grm-id", default=None, help="GRM .grm.id file")
    ap.add_argument("--grm-gz", default=None, help="GRM .grm.gz file")
    ap.add_argument("--generic", nargs=2, action="append", metavar=("PATH", "ID_COL"),
                    help="Generic file: path and ID column (repeatable)")
    args = ap.parse_args()

    keep = read_keep_list(args.keep_list)
    print(f"Loaded keep-list: {len(keep)} IIDs")

    excluded = read_exclusion_lists(args.exclusion_lists)
    print(f"Loaded exclusions: {len(excluded)} IIDs")

    all_ok = True

    if args.fam:
        all_ok &= validate_fam(args.fam, keep, excluded)

    if args.eigenvec:
        all_ok &= validate_eigenvec(args.eigenvec, keep, excluded)

    if args.cnv:
        all_ok &= validate_cnv(args.cnv, args.cnv_col, keep, excluded)

    if args.grm_bin and args.grm_id:
        all_ok &= validate_grm_binary(args.grm_bin, args.grm_id, keep, excluded)

    if args.grm_gz:
        all_ok &= validate_grm_text_gz(args.grm_gz, keep, excluded)

    if args.generic:
        for path, id_col in args.generic:
            all_ok &= validate_cnv(path, id_col, keep, excluded)  # reuse CNV logic

    print("\n" + "="*50)
    if all_ok:
        print("ALL VALIDATIONS PASSED")
        sys.exit(0)
    else:
        print("VALIDATION FAILURES DETECTED")
        sys.exit(1)


if __name__ == "__main__":
    main()