#!/usr/bin/env python3
"""Collect native CTSLEB PRS result files without selecting a best model."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


FIELDS = [
    "replicate",
    "method",
    "ancestry",
    "parameter",
    "r2",
    "adjusted_r2",
    "n",
    "status",
    "source_file",
    "note",
]


def number(value: str | None) -> str:
    """Return a finite-looking numeric field or an empty string."""
    if value is None:
        return ""
    value = value.strip()
    try:
        float(value)
    except ValueError:
        return ""
    return value


def add_row(rows: list[dict[str, str]], **values: str) -> None:
    row = {field: "" for field in FIELDS}
    row.update(values)
    rows.append(row)


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def read_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def collect_single(rows: list[dict[str, str]], replicate: str, root: Path) -> None:
    single = root / "single_ancestry_AFR"
    specs = [
        ("C+T", single / "prs_pipeline" / "CT" / "C+T_heldout_metrics.tsv"),
        ("LDpred2-inf", single / "prs_pipeline" / "LDpred2" / "LDpred2-inf_heldout_metrics.tsv"),
        ("LDpred2-grid", single / "prs_pipeline" / "LDpred2" / "LDpred2-grid_heldout_metrics.tsv"),
        ("lassosum2", single / "prs_pipeline" / "lassosum2" / "lassosum2_heldout_metrics.tsv"),
    ]
    for method, path in specs:
        if not path.is_file():
            add_row(
                rows,
                replicate=replicate,
                method=method,
                ancestry="AFR",
                status="missing_fair_metric",
                source_file=str(path),
                note="Rerun this method after the shared tuning/testing split update.",
            )
            continue
        try:
            entries = read_tsv_rows(path)
        except (OSError, csv.Error) as error:
            add_row(rows, replicate=replicate, method=method, ancestry="AFR", status="unreadable", source_file=str(path), note=str(error))
            continue
        for index, entry in enumerate(entries, start=1):
            parameter = entry.get("method", f"row_{index}")
            add_row(
                rows,
                replicate=replicate,
                method=entry.get("method", method),
                ancestry="AFR",
                parameter=parameter,
                r2=number(entry.get("primary_metric", entry.get("incremental_r2"))),
                n=number(entry.get("test_n")),
                status="complete",
                source_file=str(path),
                note="Selected using the shared tuning set and evaluated once in the shared test set.",
            )


def collect_prscsx(rows: list[dict[str, str]], replicate: str, root: Path) -> None:
    prscsx = root / "prs_pipeline" / "PRScsx"
    # AFR is the benchmark target population. EUR is the auxiliary discovery
    # population and is not a comparably held-out evaluation in this design.
    for ancestry in ("AFR",):
        path = prscsx / f"{ancestry}_PRS_sscore_metrics.tsv"
        if not path.is_file():
            add_row(rows, replicate=replicate, method="PRS-CSx", ancestry=ancestry, status="missing", source_file=str(path))
            continue
        try:
            # PRS-CSx uses the common evaluator, which writes a TSV despite
            # the historic .tsv/.csv naming mixed elsewhere in the pipeline.
            entries = read_tsv_rows(path)
        except (OSError, csv.Error) as error:
            add_row(rows, replicate=replicate, method="PRS-CSx", ancestry=ancestry, status="unreadable", source_file=str(path), note=str(error))
            continue
        for entry in entries:
            add_row(
                rows,
                replicate=replicate,
                method="PRS-CSx",
                ancestry=ancestry,
                parameter="auto",
                r2=number(entry.get("r2")),
                adjusted_r2=number(entry.get("adjusted_r2")),
                n=number(entry.get("n")),
                status="complete",
                source_file=str(path),
                note="Evaluated in the shared held-out test set.",
            )


def collect_prosper(rows: list[dict[str, str]], replicate: str, root: Path) -> None:
    path = root / "method_runs" / "multi_prosper" / "after_ensemble_AFR" / "R2.txt"
    prs_file = path.with_name("PROSPER_prs_file.txt")
    metrics = path.with_name("PROSPER_fair_metrics.tsv")
    if not path.is_file() or not prs_file.is_file() or not metrics.is_file():
        missing = []
        if not path.is_file():
            missing.append("R2.txt")
        if not prs_file.is_file():
            missing.append("PROSPER_prs_file.txt")
        if not metrics.is_file():
            missing.append("PROSPER_fair_metrics.tsv")
        add_row(rows, replicate=replicate, method="PROSPER", ancestry="AFR", status="missing", source_file=str(path), note="Missing " + ", ".join(missing))
        return
    try:
        entries = read_tsv_rows(metrics)
    except (OSError, csv.Error) as error:
        add_row(rows, replicate=replicate, method="PROSPER", ancestry="AFR", status="unreadable", source_file=str(metrics), note=str(error))
        return
    for entry in entries:
        add_row(
            rows,
            replicate=replicate,
            method=entry.get("method", "PROSPER"),
            ancestry="AFR",
            r2=number(entry.get("test_r2")),
            n=number(entry.get("test_n")),
            status="complete",
            source_file=str(metrics),
            note="Selected using the shared tuning set and evaluated once in the shared test set.",
        )


def collect_ctsleb(rows: list[dict[str, str]], replicate: str, root: Path) -> None:
    method_root = root / "method_runs" / "multi_ctsleb"
    done = root / "method_runs" / "multi_ctsleb.done"
    metrics = method_root / "AFR_EUR_ctsleb_fair_metrics.tsv"

    required = [done, metrics]
    if not all(path.is_file() for path in required):
        missing = ", ".join(path.name for path in required if not path.is_file())
        add_row(
            rows,
            replicate=replicate,
            method="CT-SLEB",
            ancestry="AFR",
            status="missing",
            source_file=str(method_root),
            note="Missing " + missing,
        )
        return

    try:
        entries = read_tsv_rows(metrics)
    except (OSError, csv.Error) as error:
        add_row(rows, replicate=replicate, method="CT-SLEB", ancestry="AFR", status="unreadable", source_file=str(metrics), note=str(error))
        return
    for entry in entries:
        add_row(
            rows,
            replicate=replicate,
            method=entry.get("method", "CT-SLEB"),
            ancestry="AFR",
            parameter=entry.get("selected_prs_column", ""),
            r2=number(entry.get("test_r2")),
            n=number(entry.get("test_n")),
            status="complete",
            source_file=str(metrics),
            note="Selected using the shared tuning set and evaluated once in the shared test set.",
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--replicates", required=True, type=int)
    args = parser.parse_args()
    if args.replicates < 1:
        parser.error("--replicates must be positive")

    results: list[dict[str, str]] = []
    for index in range(1, args.replicates + 1):
        replicate = f"rep{index:02d}"
        root = args.output_dir / "prs_inputs_ctsleb" / "AFR_EUR" / replicate
        collect_single(results, replicate, root)
        collect_prscsx(results, replicate, root)
        collect_prosper(results, replicate, root)
        collect_ctsleb(results, replicate, root)

    destination = args.output_dir / "benchmark_results_raw.tsv"
    with destination.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)
    print(f"Wrote {destination} with {len(results)} rows.")


if __name__ == "__main__":
    main()
