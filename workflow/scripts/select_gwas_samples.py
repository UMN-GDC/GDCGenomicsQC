#!/usr/bin/env python3
"""Split samples into a discovery (GWAS) and target (PRS) set.

The split is deterministic for a given seed: samples are sorted by (FID, IID)
before shuffling, so the result does not depend on input order or the Python
version. When splitting is disabled (no PRS tie-in) both keep files contain
every sample.
"""

import argparse
import random
import sys


def read_fam(path):
    ids = []
    with open(path) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) >= 2:
                ids.append((fields[0], fields[1]))
    return ids


def write_keep(path, ids):
    with open(path, "w") as handle:
        for fid, iid in ids:
            handle.write(f"{fid}\t{iid}\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fam", required=True, help="PLINK .fam file")
    parser.add_argument("--gwas-keep", required=True)
    parser.add_argument("--prs-keep", required=True)
    parser.add_argument("--fraction", type=float, default=0.5,
                        help="Fraction of samples used for GWAS discovery")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--split", action="store_true",
                        help="Create disjoint discovery/target sets")
    args = parser.parse_args()

    ids = read_fam(args.fam)
    if not ids:
        sys.exit(f"No samples found in {args.fam}")

    if not args.split:
        write_keep(args.gwas_keep, ids)
        write_keep(args.prs_keep, ids)
        print(f"Split disabled: {len(ids)} samples in both discovery and target sets")
        return

    if len(ids) < 4:
        sys.exit("At least 4 samples are required to split discovery/target sets")

    if not 0 < args.fraction < 1:
        sys.exit("--fraction must be strictly between 0 and 1 when --split is used")

    shuffled = sorted(ids)
    random.Random(args.seed).shuffle(shuffled)

    n_gwas = round(len(shuffled) * args.fraction)
    n_gwas = max(2, min(n_gwas, len(shuffled) - 2))

    write_keep(args.gwas_keep, shuffled[:n_gwas])
    write_keep(args.prs_keep, shuffled[n_gwas:])
    print(
        f"Discovery (GWAS): {n_gwas} samples; target (PRS): "
        f"{len(shuffled) - n_gwas} samples (seed={args.seed})"
    )


if __name__ == "__main__":
    main()
