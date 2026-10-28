#!/usr/bin/env python3
"""Prepare per-subset phenotype/covariate files for GWAS and PRS.

Reads a user phenotype table (``FID IID phenotype [phenotype ...]``) and an
optional quantitative covariate table, selects one phenotype and the first
``--npc`` principal components, and writes FID/IID-restricted tables for the
discovery (GWAS) and target (PRS) sample sets.
"""

import argparse
import sys


def split_line(line):
    return line.rstrip("\n").rstrip("\r").split("\t") if "\t" in line else line.split()


def read_table(path):
    with open(path) as handle:
        lines = [line for line in handle if line.strip()]
    if not lines:
        sys.exit(f"Empty table: {path}")
    header = [h.lstrip("#") for h in split_line(lines[0])]
    rows = [split_line(line) for line in lines[1:]]
    return header, rows


def id_columns(header):
    lower = [h.strip().lower() for h in header]
    iid = lower.index("iid") if "iid" in lower else 0
    fid = lower.index("fid") if "fid" in lower else None
    return fid, iid


def index_of(header, value, id_cols, label):
    ids = {v for v in id_cols if v is not None}
    traits = [i for i in range(len(header)) if i not in ids]
    if isinstance(value, int):
        if value < 1 or value > len(traits):
            sys.exit(f"{label} index {value} out of range (1..{len(traits)})")
        return traits[value - 1]
    lookup = {h.strip().lower(): i for i, h in enumerate(header)}
    key = str(value).strip().lower()
    if key not in lookup:
        sys.exit(f"{label} column '{value}' not found in header: {header}")
    return lookup[key]


def read_keep(path):
    keep = set()
    with open(path) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) >= 2:
                keep.add((fields[0], fields[1]))
            elif fields:
                keep.add((fields[0], fields[0]))
    return keep


def sample_key(fid, iid, use_fid):
    return (fid, iid) if use_fid else (iid, iid)


def write_table(path, header, rows):
    with open(path, "w") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join(row) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pheno", required=True)
    parser.add_argument("--mpheno", default="1",
                        help="Phenotype column name or 1-based trait index")
    parser.add_argument("--eigenvec", default=None,
                        help="PLINK2/GCTA eigenvec file (PC columns)")
    parser.add_argument("--npc", type=int, default=10)
    parser.add_argument("--qcovar", default=None,
                        help="Optional quantitative covariate table (FID IID ...)")
    parser.add_argument("--gwas-keep", required=True)
    parser.add_argument("--prs-keep", required=True)
    parser.add_argument("--gwas-pheno", required=True)
    parser.add_argument("--prs-pheno", required=True)
    parser.add_argument("--gwas-covar", required=True)
    parser.add_argument("--prs-covar", required=True)
    parser.add_argument("--prs-pcs", required=True)
    parser.add_argument("--tractor-pheno", required=True)
    args = parser.parse_args()

    header, rows = read_table(args.pheno)
    fid_idx, iid_idx = id_columns(header)
    use_fid = fid_idx is not None
    trait_arg = args.mpheno
    if isinstance(trait_arg, str) and trait_arg.isdigit():
        trait_arg = int(trait_arg)
    trait_idx = index_of(header, trait_arg, (fid_idx, iid_idx), "Phenotype")

    pc_names = []
    pc_rows = {}
    if args.eigenvec and args.npc > 0:
        pc_header, pc_data = read_table(args.eigenvec)
        pfid, piid = id_columns(pc_header)
        pc_cols = [i for i, name in enumerate(pc_header) if name.lower().startswith("pc")]
        pc_cols = pc_cols[: args.npc]
        if not pc_cols:
            sys.exit(f"No PC columns found in {args.eigenvec}")
        pc_names = [pc_header[i] for i in pc_cols]
        for row in pc_data:
            key = sample_key(row[pfid], row[piid], use_fid)
            pc_rows[key] = [row[i] for i in pc_cols]

    qcovar_names = []
    qcovar_rows = {}
    if args.qcovar:
        qheader, qdata = read_table(args.qcovar)
        qfid, qiid = id_columns(qheader)
        skip = {i for i in (qfid, qiid) if i is not None}
        qcols = [i for i in range(len(qheader)) if i not in skip]
        qcovar_names = [qheader[i] for i in qcols]
        for row in qdata:
            key = sample_key(row[qfid], row[qiid], use_fid)
            qcovar_rows[key] = [row[i] for i in qcols]

    gwas_keep = read_keep(args.gwas_keep)
    prs_keep = read_keep(args.prs_keep)
    covar_names = pc_names + qcovar_names

    def select(keep, want_pheno):
        out = []
        for row in rows:
            fid = row[fid_idx] if use_fid else row[iid_idx]
            iid = row[iid_idx]
            key = sample_key(fid, iid, use_fid)
            if key not in keep:
                continue
            pheno = row[trait_idx] if trait_idx < len(row) else "NA"
            if pheno in ("", "NA", "-9", "nan", "NaN"):
                continue
            covars = []
            missing = False
            for name in pc_names:
                value = pc_rows.get(key, [None] * len(pc_names))[pc_names.index(name)]
                if value in (None, "", "NA"):
                    missing = True
                    break
                covars.append(value)
            if not missing:
                for name in qcovar_names:
                    value = qcovar_rows.get(key, [None] * len(qcovar_names))[
                        qcovar_names.index(name)
                    ]
                    if value in (None, "", "NA"):
                        missing = True
                        break
                    covars.append(value)
            if missing:
                continue
            out.append([fid, iid, pheno] + covars)
        return out

    gwas_rows = select(gwas_keep, True)
    prs_rows = select(prs_keep, False)
    if not gwas_rows:
        sys.exit("No discovery samples remain after intersecting phenotype/covariates")

    write_table(args.gwas_pheno, ["FID", "IID", "PHENO"], [r[:3] for r in gwas_rows])
    write_table(args.tractor_pheno, ["FID", "IID", "PHENO"] + covar_names,
                [r[:3] + r[3:] for r in gwas_rows])
    write_table(args.gwas_covar, ["FID", "IID"] + covar_names,
                [r[:2] + r[3:] for r in gwas_rows])

    if prs_rows:
        write_table(args.prs_pheno, ["FID", "IID", "PHENO"], [r[:3] for r in prs_rows])
        write_table(args.prs_covar, ["FID", "IID"] + covar_names,
                    [r[:2] + r[3:] for r in prs_rows])
    else:
        write_table(args.prs_pheno, ["FID", "IID", "PHENO"], [])
        write_table(args.prs_covar, ["FID", "IID"] + covar_names, [])

    write_table(args.prs_pcs, ["FID", "IID"] + pc_names,
                [r[:2] + r[3:3 + len(pc_names)] for r in prs_rows])

    with open(args.gwas_covar + ".columns", "w") as handle:
        handle.write("\n".join(covar_names) + ("\n" if covar_names else ""))

    print(f"Discovery samples: {len(gwas_rows)}; target samples: {len(prs_rows)}")
    print(f"Covariates ({len(covar_names)}): {', '.join(covar_names) or 'none'}")


if __name__ == "__main__":
    main()
