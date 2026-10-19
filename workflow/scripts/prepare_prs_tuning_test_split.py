#!/usr/bin/env python3
"""Create an auditable, optionally family-preserving PRS tuning/test split.

The input phenotype table must contain FID, IID, and PHENO.  A group file is
optional, but should be supplied for cohorts containing related participants:
members with the same group value are always assigned to the same partition.
"""

from __future__ import annotations

import argparse
import csv
import json
import random
from collections import defaultdict
from pathlib import Path


def read_tsv(path: Path, required: set[str]) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = set(reader.fieldnames or [])
        missing = sorted(required - columns)
        if missing:
            raise ValueError(f"{path} is missing required column(s): {', '.join(missing)}")
        return list(reader)


def identifier(row: dict[str, str]) -> tuple[str, str]:
    return str(row["FID"]), str(row["IID"])


def write_table(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_keep(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for row in rows:
            writer.writerow([row["FID"], row["IID"]])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a reproducible, optionally family-preserving PRS tuning/test split."
    )
    parser.add_argument("--pheno", required=True, type=Path, help="Headered FID/IID/PHENO TSV.")
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--seed", required=True, type=int)
    parser.add_argument("--tuning-fraction", type=float, default=0.5)
    parser.add_argument("--label", default="target")
    parser.add_argument(
        "--group-file",
        type=Path,
        help="Optional headered FID/IID/GROUP TSV. Members of a GROUP stay together.",
    )
    parser.add_argument("--group-column", default="GROUP")
    parser.add_argument(
        "--covariates",
        type=Path,
        help="Optional headered FID/IID covariate TSV to partition alongside the phenotype.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not 0 < args.tuning_fraction < 1:
        raise SystemExit("--tuning-fraction must be strictly between 0 and 1.")

    phenotype_rows = read_tsv(args.pheno, {"FID", "IID", "PHENO"})
    if not phenotype_rows:
        raise SystemExit("Phenotype table is empty.")
    if len({identifier(row) for row in phenotype_rows}) != len(phenotype_rows):
        raise SystemExit("Phenotype table contains duplicate FID/IID pairs.")

    group_by_id: dict[tuple[str, str], str] = {}
    group_source = "individual"
    if args.group_file:
        group_rows = read_tsv(args.group_file, {"FID", "IID", args.group_column})
        group_by_id = {
            identifier(row): str(row[args.group_column])
            for row in group_rows
            if str(row[args.group_column]).strip()
        }
        missing_group = [identifier(row) for row in phenotype_rows if identifier(row) not in group_by_id]
        if missing_group:
            raise SystemExit(
                f"{len(missing_group)} phenotype participants are missing a family/group assignment."
            )
        group_source = str(args.group_file)

    groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in phenotype_rows:
        fid, iid = identifier(row)
        group = group_by_id.get((fid, iid), f"individual::{fid}::{iid}")
        groups[group].append(row)

    # Shuffle groups, never individual rows, so related participants cannot leak.
    group_ids = list(groups)
    random.Random(args.seed).shuffle(group_ids)
    target_n = round(len(phenotype_rows) * args.tuning_fraction)
    tuning_groups: set[str] = set()
    tuning_n = 0
    for group in group_ids:
        if tuning_n == 0 or abs((tuning_n + len(groups[group])) - target_n) <= abs(tuning_n - target_n):
            tuning_groups.add(group)
            tuning_n += len(groups[group])

    tuning_rows = [row for group in group_ids if group in tuning_groups for row in groups[group]]
    testing_rows = [row for group in group_ids if group not in tuning_groups for row in groups[group]]
    if not tuning_rows or not testing_rows:
        raise SystemExit("Split produced an empty tuning or testing cohort; revise the tuning fraction.")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    columns = ["FID", "IID", "PHENO"]
    tuning_pheno = args.out_dir / f"{args.label}_tuning.pheno.tsv"
    testing_pheno = args.out_dir / f"{args.label}_testing.pheno.tsv"
    tuning_keep = args.out_dir / f"{args.label}_tuning.keep.txt"
    testing_keep = args.out_dir / f"{args.label}_testing.keep.txt"
    write_table(tuning_pheno, tuning_rows, columns)
    write_table(testing_pheno, testing_rows, columns)
    write_keep(tuning_keep, tuning_rows)
    write_keep(testing_keep, testing_rows)

    covariate_counts: dict[str, int] = {}
    if args.covariates:
        covariate_rows = read_tsv(args.covariates, {"FID", "IID"})
        covariate_by_id = {identifier(row): row for row in covariate_rows}
        missing_covariates = [identifier(row) for row in phenotype_rows if identifier(row) not in covariate_by_id]
        if missing_covariates:
            raise SystemExit(
                f"{len(missing_covariates)} phenotype participants are missing covariates; "
                "create a complete analysis table before splitting."
            )
        covariate_columns = list(covariate_rows[0])
        tuning_covariates = [covariate_by_id[identifier(row)] for row in tuning_rows]
        testing_covariates = [covariate_by_id[identifier(row)] for row in testing_rows]
        write_table(args.out_dir / f"{args.label}_tuning.covariates.tsv", tuning_covariates, covariate_columns)
        write_table(args.out_dir / f"{args.label}_testing.covariates.tsv", testing_covariates, covariate_columns)
        covariate_counts = {"tuning_covariates": len(tuning_covariates), "testing_covariates": len(testing_covariates)}

    assignments = []
    for group in group_ids:
        partition = "tuning" if group in tuning_groups else "testing"
        for row in groups[group]:
            assignments.append({"FID": row["FID"], "IID": row["IID"], "GROUP": group, "PARTITION": partition})
    write_table(args.out_dir / f"{args.label}_split_assignments.tsv", assignments, ["FID", "IID", "GROUP", "PARTITION"])

    manifest = {
        "label": args.label,
        "phenotype_file": str(args.pheno),
        "group_source": group_source,
        "seed": args.seed,
        "requested_tuning_fraction": args.tuning_fraction,
        "n_participants": len(phenotype_rows),
        "n_groups": len(groups),
        "n_tuning": len(tuning_rows),
        "n_testing": len(testing_rows),
        "n_tuning_groups": len(tuning_groups),
        "n_testing_groups": len(groups) - len(tuning_groups),
        "files": {
            "tuning_keep": str(tuning_keep),
            "testing_keep": str(testing_keep),
            "tuning_pheno": str(tuning_pheno),
            "testing_pheno": str(testing_pheno),
        },
        **covariate_counts,
    }
    (args.out_dir / f"{args.label}_split_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
