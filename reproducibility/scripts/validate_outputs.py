#!/usr/bin/env python3
"""Validate the public GDCGenomicsQC v1 QC smoke-test output contract.

This checks output existence and exact sample/variant metadata against the
public synthetic fixture. It does *not* verify binary .pgen genotype calls,
QC threshold exclusion behavior, or optional downstream modules.

Backward compatible with run_smoke_test.sh:
  python validate_outputs.py --output-dir DIR --expectations EXPECTED_JSON
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


class ValidationError(Exception):
    """An expected smoke-test output or invariant was not satisfied."""


def load_json(path: Path) -> dict:
    try:
        obj = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValidationError(f"Cannot read JSON {path}: {exc}") from exc
    if not isinstance(obj, dict):
        raise ValidationError(f"Expected a JSON object: {path}")
    return obj


def parse_plink_table(path: Path, required_columns: tuple[str, ...]) -> list[dict[str, str]]:
    """Read PLINK2 whitespace-delimited .psam/.pvar, including # headers.

    The parser skips only '##' metadata (not '#IID' or '#CHROM' headers).
    """
    header: list[str] | None = None
    records: list[dict[str, str]] = []
    try:
        with path.open(encoding="utf-8") as handle:
            for number, raw in enumerate(handle, start=1):
                line = raw.strip()
                if not line or line.startswith("##"):
                    continue
                fields = line.split()
                if header is None:
                    if not fields[0].startswith("#"):
                        raise ValidationError(f"{path}:{number}: missing PLINK header")
                    header = [name.lstrip("#") for name in fields]
                    if len(header) != len(set(header)):
                        raise ValidationError(f"{path}:{number}: duplicate header columns")
                    missing = set(required_columns) - set(header)
                    if missing:
                        raise ValidationError(f"{path}:{number}: missing columns {sorted(missing)}")
                    continue
                if line.startswith("#"):
                    raise ValidationError(f"{path}:{number}: unexpected comment/header in records")
                if len(fields) != len(header):
                    raise ValidationError(
                        f"{path}:{number}: expected {len(header)} columns, found {len(fields)}"
                    )
                records.append(dict(zip(header, fields)))
    except (OSError, UnicodeError) as exc:
        raise ValidationError(f"Cannot read {path}: {exc}") from exc
    if header is None:
        raise ValidationError(f"{path}: no PLINK header found")
    return records


def assert_exact_ids(observed: list[str], expected: list[str], name: str) -> None:
    if len(observed) != len(set(observed)):
        raise ValidationError(f"{name}: duplicate identifiers detected")
    if len(expected) != len(set(expected)):
        raise ValidationError(f"{name}: duplicate identifiers in fixture truth")
    observed_set, expected_set = set(observed), set(expected)
    if observed_set != expected_set:
        missing = sorted(expected_set - observed_set)
        unexpected = sorted(observed_set - expected_set)
        raise ValidationError(
            f"{name}: expected {len(expected_set)} identifiers, observed {len(observed_set)}; "
            f"missing={missing[:10]} unexpected={unexpected[:10]}"
        )
    print(f"  PASS {name}: exact {len(expected)} identifiers (no duplicates)")


def validate(output_dir: Path, expectations: Path, truth_path: Path) -> None:
    expected = load_json(expectations)
    truth = load_json(truth_path)

    if expected.get("fixture_version") != truth.get("fixture_version"):
        raise ValidationError("Expected-output contract and truth file use different fixture versions")
    if expected.get("target") != "run_initialQC":
        raise ValidationError("This validator is intended for the run_initialQC smoke-test target")

    required_files = expected.get("required_files")
    if not isinstance(required_files, list) or not required_files:
        raise ValidationError("expected_outputs.json must contain a nonempty required_files list")
    for name in required_files:
        if not isinstance(name, str) or Path(name).is_absolute() or ".." in Path(name).parts:
            raise ValidationError(f"Unsafe or invalid required output path: {name!r}")
    missing = [name for name in required_files if not (output_dir / name).is_file()]
    if missing:
        raise ValidationError("Missing expected files:\n  " + "\n  ".join(missing))
    empty = [name for name in required_files if (output_dir / name).stat().st_size == 0]
    if empty:
        raise ValidationError("Empty expected files:\n  " + "\n  ".join(empty))
    print(f"  PASS required files: all {len(required_files)} present and nonempty")

    samples = parse_plink_table(output_dir / "full" / "f1.psam", ("IID",))
    variants = parse_plink_table(
        output_dir / "full" / "f1.pvar", ("CHROM", "POS", "ID", "REF", "ALT")
    )
    truth_samples = truth.get("samples")
    truth_variants = truth.get("expected_variant_ids")
    if not isinstance(truth_samples, list) or not all(isinstance(x, str) for x in truth_samples):
        raise ValidationError("truth.json must contain samples as a list of strings")
    if not isinstance(truth_variants, list) or not all(isinstance(x, str) for x in truth_variants):
        raise ValidationError("truth.json must contain expected_variant_ids as a list of strings")

    limits = expected.get("minimum_counts", {})
    if len(samples) < limits.get("samples", 0) or len(variants) < limits.get("variants", 0):
        raise ValidationError(
            f"Below minimum counts: samples={len(samples)}, variants={len(variants)}"
        )
    if len(samples) != truth.get("sample_count") or len(variants) != truth.get("variant_count"):
        raise ValidationError(
            f"Fixture counts differ: observed {len(samples)} samples/{len(variants)} variants; "
            f"truth {truth.get('sample_count')}/{truth.get('variant_count')}"
        )

    assert_exact_ids([row["IID"] for row in samples], truth_samples, "sample IDs")
    assert_exact_ids([row["ID"] for row in variants], truth_variants, "variant IDs")

    # The public fixture uses chr:pos:ref:alt IDs. Verify that each ID also
    # agrees with the corresponding metadata columns, not just with the truth list.
    for record in variants:
        # Normalize chromosome labels
        chrom = record["CHROM"].removeprefix("chr")
        variant_id = record["ID"]

        parts = variant_id.split(":")
        if len(parts) != 4:
            raise ValidationError(
                f"Invalid variant ID format: {variant_id!r}"
            )

        id_chrom = parts[0].removeprefix("chr")

        reconstructed = ":".join(
            [chrom, record["POS"], record["REF"], record["ALT"]]
        )
        normalized_id = ":".join([id_chrom] + parts[1:])

        if normalized_id != reconstructed:
            raise ValidationError(
                f"Variant metadata mismatch: ID={variant_id!r}, "
                f"CHROM:POS:REF:ALT={reconstructed!r}"
            )


    print("  PASS variant metadata: IDs match CHROM:POS:REF:ALT")
    print(f"PASS: QC smoke-test output contract ({len(samples)} samples, {len(variants)} variants)")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--expectations", type=Path, required=True)
    parser.add_argument(
        "--truth", type=Path, default=None,
        help="Fixture truth JSON (default: fixtures/truth.json next to expected_outputs.json)",
    )
    args = parser.parse_args()
    truth_path = args.truth or (args.expectations.parent / "fixtures" / "truth.json")
    try:
        validate(args.output_dir, args.expectations, truth_path)
    except (ValidationError, OSError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
