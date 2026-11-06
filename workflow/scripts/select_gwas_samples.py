#!/usr/bin/env python3
"""Write the complete sample set from a PLINK .fam file for GWAS."""

import argparse
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
    args = parser.parse_args()

    ids = read_fam(args.fam)
    if not ids:
        sys.exit(f"No samples found in {args.fam}")

    write_keep(args.gwas_keep, ids)
    print(f"Selected {len(ids)} samples for GWAS")


if __name__ == "__main__":
    main()
