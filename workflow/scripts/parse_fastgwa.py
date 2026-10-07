#!/usr/bin/env python3
"""Parse h2 (Heritability), Vg, Ve, Vp and p-value from a fastGWA log file."""
import argparse
import re
import sys
from pathlib import Path

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--log", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    text = Path(args.log).read_text(errors="replace")
    m = re.search(r"Heritability\s*=\s*([0-9eE+\-.]+)\s*\(Pval\s*=\s*([0-9eE+\-.]+)\)", text)
    if not m:
        sys.exit(f"could not find Heritability line in {args.log}")
    h2, pval = float(m.group(1)), float(m.group(2))
    comp = {}
    for name in ("Vg", "Ve", "Vp"):
        mm = re.search(rf"^{name}\s+([0-9eE+\-.]+)\s+([0-9eE+\-.]+)?", text, re.M)
        if mm:
            comp[f"{name}"] = float(mm.group(1))
            if mm.group(2):
                comp[f"{name}_se"] = float(mm.group(2))
    row = {"h2": h2, "pval": pval, **comp}
    pd.DataFrame([row]).to_csv(args.out, index=False)


if __name__ == "__main__":
    main()
