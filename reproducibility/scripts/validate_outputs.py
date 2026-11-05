#!/usr/bin/env python3
"""Validate the portable smoke-test outputs without scientific overclaiming."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def rows(path: Path) -> int:
    return sum(1 for line in path.read_text(encoding="utf-8").splitlines()
               if line and not line.startswith("#"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--expectations", type=Path, required=True)
    args = parser.parse_args()
    expected = json.loads(args.expectations.read_text(encoding="utf-8"))
    missing = [name for name in expected["required_files"] if not (args.output_dir / name).is_file()]
    if missing:
        sys.exit("Missing expected files:\n  " + "\n  ".join(missing))
    samples = rows(args.output_dir / "full" / "f1.psam")
    variants = rows(args.output_dir / "full" / "f1.pvar")
    limits = expected["minimum_counts"]
    if samples < limits["samples"] or variants < limits["variants"]:
        sys.exit(f"QC retained too little data: samples={samples}, variants={variants}")
    print(f"PASS: {samples} samples and {variants} variants passed the portable smoke test.")


if __name__ == "__main__":
    main()
