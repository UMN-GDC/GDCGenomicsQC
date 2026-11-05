#!/usr/bin/env python3
"""Evaluate v1 global-ancestry classifications against a public truth file."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path


def read_truth(path: Path) -> dict[str, str]:
    truth: dict[str, str] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()
            if not fields or fields[0].lower() in {"fid", "familyid"}:
                continue
            if len(fields) < 3:
                raise ValueError(f"Expected FID, IID, and ancestry in {path}: {line!r}")
            truth[fields[1]] = fields[-1]
    if not truth:
        raise ValueError(f"No truth labels found in {path}")
    return truth


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--classifications", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--model", default="pca")
    args = parser.parse_args()

    truth = read_truth(args.truth)
    prediction_col = f"{args.model}_predicted"
    confidence_col = f"{args.model}_confidence"
    observed: list[tuple[str, str, str, float | None]] = []
    with args.classifications.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"IID", prediction_col, confidence_col}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"Expected columns {sorted(required)} in {args.classifications}")
        for row in reader:
            iid = row["IID"]
            if iid not in truth:
                continue
            confidence = float(row[confidence_col]) if row[confidence_col] else None
            observed.append((iid, truth[iid], row[prediction_col], confidence))

    if not observed:
        raise ValueError("No classified IIDs matched the truth file")

    labels = sorted(set(truth.values()) | {predicted for _, _, predicted, _ in observed})
    matrix: dict[str, Counter[str]] = defaultdict(Counter)
    assigned = [row for row in observed if row[2] != "uncertain"]
    correct_all = sum(actual == predicted for _, actual, predicted, _ in observed)
    correct_assigned = sum(actual == predicted for _, actual, predicted, _ in assigned)
    for _, actual, predicted, _ in observed:
        matrix[actual][predicted] += 1

    args.output_dir.mkdir(parents=True, exist_ok=True)
    with (args.output_dir / "confusion_matrix.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["truth\\prediction", *labels])
        for actual in sorted(matrix):
            writer.writerow([actual, *(matrix[actual][predicted] for predicted in labels)])

    summary = {
        "model": args.model,
        "truth_samples": len(truth),
        "matched_samples": len(observed),
        "assigned_samples": len(assigned),
        "uncertain_samples": len(observed) - len(assigned),
        "uncertain_rate": (len(observed) - len(assigned)) / len(observed),
        "accuracy_all_samples": correct_all / len(observed),
        "accuracy_assigned_samples": correct_assigned / len(assigned) if assigned else None,
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
