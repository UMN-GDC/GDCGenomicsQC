#!/usr/bin/env python3
"""Convert GWAS output (PLINK2, GCTA fastGWA, TRACTOR) to PRS-ready sumstats.

Output columns (tab-separated, one header line)::

    SNP  CHR  POS  A1  A2  BETA  SE  P  N

``A1`` is the effect allele. Duplicate variant IDs are dropped, keeping the
first occurrence. Tractor output is selected with ``--ancestry``.
"""

import argparse
import sys


def rows(path):
    with open(path) as handle:
        header = None
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                continue
            if header is None:
                header = line.rstrip("\n").split("\t")
                continue
            yield header, line.rstrip("\n").split("\t")


def column(header, *names):
    lookup = {name.strip(): i for i, name in enumerate(header)}
    for name in names:
        if name in lookup:
            return lookup[name]
    sys.exit(f"None of {names} found in header: {header}")


def clean(value):
    return value if value not in ("", "NA", "nan", "NaN", ".") else None


def convert_plink2(path, out):
    emitted = 0
    seen = set()
    with open(out, "w") as handle:
        handle.write("SNP\tCHR\tPOS\tA1\tA2\tBETA\tSE\tP\tN\n")
        for header, row in rows(path):
            idx = {name: i for i, name in enumerate(header)}
            if idx.get("TEST") is not None and row[idx["TEST"]] != "ADD":
                continue
            snp = clean(row[idx["ID"]])
            chr_ = clean(row[idx["CHROM"]])
            pos = clean(row[idx["POS"]])
            ref = clean(row[idx["REF"]])
            alt = clean(row[idx["ALT"]])
            a1 = clean(row[idx["A1"]])
            beta = clean(row[idx["BETA"]])
            se = clean(row[idx["SE"]])
            p = clean(row[idx["P"]])
            n = clean(row[idx["OBS_CT"]])
            if None in (snp, chr_, a1, alt, ref, beta, se, p, n) or snp in seen:
                continue
            a2 = ref if a1 == alt else alt
            seen.add(snp)
            handle.write(
                f"{snp}\t{chr_}\t{pos or '0'}\t{a1}\t{a2}\t{beta}\t{se}\t{p}\t{n}\n"
            )
            emitted += 1
    return emitted


def convert_gcta(path, out):
    emitted = 0
    seen = set()
    with open(out, "w") as handle:
        handle.write("SNP\tCHR\tPOS\tA1\tA2\tBETA\tSE\tP\tN\n")
        for header, row in rows(path):
            vals = dict(zip(header, row))
            snp = clean(vals.get("SNP"))
            chr_ = clean(vals.get("CHR"))
            pos = clean(vals.get("POS"))
            a1 = clean(vals.get("A1"))
            a2 = clean(vals.get("A2"))
            beta = clean(vals.get("BETA"))
            se = clean(vals.get("SE"))
            p = clean(vals.get("P"))
            n = clean(vals.get("N"))
            if None in (snp, chr_, a1, a2, beta, se, p, n) or snp in seen:
                continue
            seen.add(snp)
            handle.write(
                f"{snp}\t{chr_}\t{pos or '0'}\t{a1}\t{a2}\t{beta}\t{se}\t{p}\t{n}\n"
            )
            emitted += 1
    return emitted


def convert_tractor(path, out, ancestry):
    suffix = f"anc{ancestry}"
    emitted = 0
    seen = set()
    with open(out, "w") as handle:
        handle.write("SNP\tCHR\tPOS\tA1\tA2\tBETA\tSE\tP\tN\n")
        for header, row in rows(path):
            vals = dict(zip(header, row))
            snp = clean(vals.get("ID"))
            chr_ = clean(vals.get("CHR"))
            pos = clean(vals.get("POS"))
            ref = clean(vals.get("REF"))
            alt = clean(vals.get("ALT"))
            beta = clean(vals.get(f"beta_{suffix}"))
            se = clean(vals.get(f"se_{suffix}"))
            p = clean(vals.get(f"pval_{suffix}"))
            n = clean(vals.get("N"))
            if None in (snp, chr_, ref, alt, beta, se, p, n) or snp in seen:
                continue
            seen.add(snp)
            handle.write(
                f"{snp}\t{chr_}\t{pos or '0'}\t{alt}\t{ref}\t{beta}\t{se}\t{p}\t{n}\n"
            )
            emitted += 1
    return emitted


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--method", required=True, choices=["plink2", "gcta", "tractor"])
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--ancestry", type=int, default=None,
                        help="Ancestry index for --method tractor")
    args = parser.parse_args()

    if args.method == "plink2":
        emitted = convert_plink2(args.input, args.output)
    elif args.method == "gcta":
        emitted = convert_gcta(args.input, args.output)
    else:
        if args.ancestry is None:
            sys.exit("--ancestry is required for --method tractor")
        emitted = convert_tractor(args.input, args.output, args.ancestry)

    if emitted == 0:
        sys.exit(f"No usable variants converted from {args.input}")
    print(f"Wrote {emitted} variants to {args.output}")


if __name__ == "__main__":
    main()
