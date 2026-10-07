#!/usr/bin/env python3
"""
Build release keep-list from exclusion sources and QC-passing subjects.

Inputs:
- --source-fam: PLINK .fam (or .psam converted) with all subjects post-QC (raw IDs)
- --identifiers: Crosswalk file mapping raw_id -> release_id (+ optional suffix)
- --exclusion-sources: One or more exclusion files (CSV/TSV/Excel)
- --par-visit: Optional file listing subjects eligible for release (raw IDs)
- --id-col: Column name for subject ID in source files
- --release-id-col: Column name for release ID in crosswalk
- --suffix-col: Optional column for C/M suffix in crosswalk

Outputs:
- --keep-list: One-column file of release IIDs (release_id + suffix) for plink2 --keep
- --temp-fam: Modified .fam with de-identified FID/IID (all QC-passing subjects, preserves row order)
- --batch-info: Optional batch info file (IID, visit, plate_number)
- --removed-individuals: List of excluded IIDs for documentation

Logic:
1. Start with all subjects in source-fam (post-QC)
2. Apply par-visit filter (if provided) - keep only eligible
3. Apply exclusion sources - remove any matching subjects
4. Map remaining raw IDs to release IDs via crosswalk
5. Output keep_list (release_id + suffix) and temp.fam (de-IDed FID/IID, all QC-passing)
"""
import argparse, sys, pandas as pd, numpy as np
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
    # default: whitespace-delimited (PLINK format)
    return pd.read_csv(p, sep=r"\s+", **kwargs)


def get_id_column(df, preferred):
    """Find ID column - try preferred, then common names."""
    for c in [preferred, "IID", "FID", "sample_id", "subject_id", "ID", "id"]:
        if c in df.columns:
            return c
    return df.columns[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-fam", required=True, help="PLINK .fam or .psam with all post-QC subjects")
    ap.add_argument("--identifiers", required=True, help="Crosswalk: raw_id -> release_id (+ suffix)")
    ap.add_argument("--exclusion-sources", nargs="*", default=[], help="Exclusion files (CSV/TSV/Excel)")
    ap.add_argument("--par-visit", default=None, help="Optional: eligible subjects (raw IDs)")
    ap.add_argument("--id-col", default="IID", help="ID column name in source files")
    ap.add_argument("--release-id-col", default="release_candid", help="Release ID column in crosswalk")
    ap.add_argument("--suffix-col", default=None, help="Optional C/M suffix column in crosswalk")
    ap.add_argument("--keep-list", required=True, help="Output: one-column release IIDs for plink2 --keep")
    ap.add_argument("--temp-fam", required=True, help="Output: de-identified .fam (all QC-passing, preserves row order)")
    ap.add_argument("--batch-info", default=None, help="Optional output: batch info (IID visit plate)")
    ap.add_argument("--removed-individuals", default=None, help="Optional output: excluded IIDs")
    args = ap.parse_args()

    # 1. Read source FAM/PSAM - all post-QC subjects
    source = read_table(args.source_fam, dtype=str)
    # Standardize column names
    if args.id_col not in source.columns:
        # Try common alternatives
        for alt in ["IID", "FID", "sample_id", "subject_id", "ID"]:
            if alt in source.columns:
                source = source.rename(columns={alt: args.id_col})
                break
    if args.id_col not in source.columns:
        sys.exit(f"ID column '{args.id_col}' not found in {args.source_fam}. Columns: {list(source.columns)}")

    all_raw_ids = source[args.id_col].astype(str).tolist()
    print(f"Loaded {len(all_raw_ids)} subjects from source fam")

    # 2. Optional par-visit filter
    if args.par_visit:
        pv = read_table(args.par_visit, dtype=str)
        pv_id = get_id_column(pv, args.id_col)
        pv_ids = set(pv[pv_id].astype(str))
        all_raw_ids = [rid for rid in all_raw_ids if rid in pv_ids]
        print(f"After par-visit filter: {len(all_raw_ids)} subjects")

    # 3. Load identifiers crosswalk
    idmap = read_table(args.identifiers, dtype=str)
    # Find raw ID column in crosswalk
    raw_col = None
    for c in ["pscid", "raw_id", "raw_ID", "source_id", "Study_ID", "ID"]:
        if c in idmap.columns:
            raw_col = c
            break
    if raw_col is None:
        sys.exit(f"Could not find raw ID column in crosswalk. Columns: {list(idmap.columns)}")

    # Build mapping: raw_id -> (release_id, suffix)
    idmap[raw_col] = idmap[raw_col].astype(str)
    idmap[args.release_id_col] = idmap[args.release_id_col].astype(str)

    if args.suffix_col and args.suffix_col in idmap.columns:
        idmap[args.suffix_col] = idmap[args.suffix_col].astype(str)
        mapping = dict(zip(idmap[raw_col], zip(idmap[args.release_id_col], idmap[args.suffix_col])))
    else:
        # Default suffix 'C' if not provided
        mapping = {row[raw_col]: (row[args.release_id_col], "C") for _, row in idmap.iterrows()}

    # 4. Apply exclusion sources
    excluded_raw = set()
    for excl_file in args.exclusion_sources:
        excl = read_table(excl_file, dtype=str)
        # Try to find ID column
        excl_id_col = get_id_column(excl, args.id_col)
        excl_ids = set(excl[excl_id_col].astype(str))
        excluded_raw.update(excl_ids)
        print(f"  Exclusion {excl_file}: {len(excl_ids)} IDs")

    print(f"Total excluded (raw): {len(excluded_raw)}")

    # 4. Filter: keep only raw IDs that are in all_raw_ids AND not excluded
    # Also must be present in crosswalk
    mappable_raw = set(mapping.keys())
    final_raw = [rid for rid in all_raw_ids if rid in mappable_raw and rid not in excluded_raw]
    excluded_final = [rid for rid in all_raw_ids if rid not in mappable_raw or rid in excluded_raw]

    print(f"Final release subjects (raw): {len(final_raw)}")
    print(f"Removed (excluded + unmappable): {len(excluded_final)}")

    # 5. Build de-identified IDs
    final_deid = [mapping[rid][0] + mapping[rid][1] for rid in final_raw]
    excluded_deid = [mapping[rid][0] + mapping[rid][1] if rid in mapping else rid for rid in excluded_final]

    # 6. Write keep-list (one column, release IIDs for plink2 --keep)
    Path(args.keep_list).parent.mkdir(parents=True, exist_ok=True)
    pd.Series(final_deid).to_csv(args.keep_list, index=False, header=False)
    print(f"Wrote keep-list: {args.keep_list} ({len(final_deid)} IIDs)")

    # 7. Build temp.fam: de-identified FID/IID for ALL QC-passing subjects (preserves row order)
    # Source FAM may have FID column; if not, use IID as FID
    if "FID" not in source.columns:
        source["FID"] = source[args.id_col]

    def deid_row(row):
        raw = str(row[args.id_col])
        if raw in mapping:
            rid, suf = mapping[raw]
            return pd.Series({"FID": rid, "IID": rid + suf})
        else:
            # Unmappable - use raw ID (will be filtered out by --keep)
            return pd.Series({"FID": raw, "IID": raw + "C"})

    temp_fam = source.copy()
    deid_cols = temp_fam.apply(deid_row, axis=1)
    temp_fam["FID"] = deid_cols["FID"]
    temp_fam["IID"] = deid_cols["IID"]

    Path(args.temp_fam).parent.mkdir(parents=True, exist_ok=True)
    # PLINK .fam format: FID IID PAT MAT SEX PHENO
    # Preserve original PAT/MAT/SEX/PHENO if present
    out_cols = ["FID", "IID"]
    for c in ["PAT", "MAT", "SEX", "PHENO"]:
        if c in temp_fam.columns:
            out_cols.append(c)
        else:
            temp_fam[c] = 0 if c != "SEX" else 1
            out_cols.append(c)
    temp_fam[out_cols].to_csv(args.temp_fam, sep=" ", index=False, header=False)
    print(f"Wrote temp.fam: {args.temp_fam} ({len(temp_fam)} rows)")

    # 8. Optional batch info
    if args.batch_info:
        batch_df = pd.DataFrame({"IID": final_deid, "visit": 1, "plate_number": 1})
        Path(args.batch_info).parent.mkdir(parents=True, exist_ok=True)
        batch_df.to_csv(args.batch_info, sep="\t", index=False)
        print(f"Wrote batch info: {args.batch_info}")

    # 9. Optional removed individuals
    if args.removed_individuals:
        Path(args.removed_individuals).parent.mkdir(parents=True, exist_ok=True)
        pd.Series(excluded_deid).to_csv(args.removed_individuals, index=False, header=False)
        print(f"Wrote removed individuals: {args.removed_individuals} ({len(excluded_deid)} IIDs)")


if __name__ == "__main__":
    main()