#!/usr/bin/env python3
"""Create a tiny, deterministic, non-human genomic test fixture.

The fixture is deliberately too small for scientific inference.  It tests the
input contract, identifier handling, basic QC, and reporting without exposing
participant-level data.  The truth file is a public oracle for test assertions.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


SAMPLES = [f"SYN{i:03d}" for i in range(1, 13)]
VARIANTS = [
    ("chr1", 101, "A", "G", ["0/0", "0/1", "0/0", "1/1", "0/1", "0/0", "0/1", "0/0", "0/1", "0/0", "0/1", "0/0"]),
    ("chr1", 205, "C", "T", ["0/1", "0/0", "0/1", "0/0", "1/1", "0/1", "0/0", "0/1", "0/0", "0/1", "0/0", "0/1"]),
    ("chr1", 318, "G", "A", ["0/0", "0/0", "0/1", "0/1", "0/0", "1/1", "0/1", "0/0", "0/1", "0/0", "0/0", "0/1"]),
    ("chr1", 422, "T", "C", ["0/1", "0/1", "0/0", "0/1", "0/0", "0/1", "1/1", "0/0", "0/1", "0/1", "0/0", "0/1"]),
    ("chr1", 537, "A", "C", ["0/0", "0/1", "0/0", "0/1", "0/1", "0/0", "0/1", "0/0", "1/1", "0/1", "0/0", "0/1"]),
    ("chr1", 641, "C", "G", ["0/1", "0/0", "0/1", "0/0", "0/1", "0/1", "0/0", "1/1", "0/1", "0/0", "0/1", "0/0"]),
]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    out = args.output_dir
    reference = out / "reference" / "1000G_highcoverage"
    reference.mkdir(parents=True, exist_ok=True)

    fasta = out / "reference" / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    sequence = list("A" * 700)
    for _, pos, ref, _, _ in VARIANTS:
        sequence[pos - 1] = ref
    fasta.write_text(">chr1\n" + "".join(sequence) + "\n", encoding="utf-8")

    vcf = out / "synthetic.vcf"
    header = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr1,length=700>",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(SAMPLES),
    ]
    rows = ["\t".join([chrom, str(pos), f"chr1:{pos}:{ref}:{alt}", ref, alt, "100", "PASS", ".", "GT", *genotypes])
            for chrom, pos, ref, alt, genotypes in VARIANTS]
    vcf.write_text("\n".join(header + rows) + "\n", encoding="utf-8")

    pvar = reference / "1000G_highCoveragephased.pvar"
    pvar.write_text("#CHROM\tPOS\tID\tREF\tALT\n" + "\n".join(
        f"{chrom}\t{pos}\tchr1:{pos}:{ref}:{alt}\t{ref}\t{alt}"
        for chrom, pos, ref, alt, _ in VARIANTS
    ) + "\n", encoding="utf-8")

    keep = out / "provided_ancestry.tsv"
    keep.write_text("IID\tancestry\n" + "\n".join(f"{sample}\tEUR" for sample in SAMPLES) + "\n", encoding="utf-8")
    truth = {
        "fixture_version": "1.0.0",
        "samples": SAMPLES,
        "sample_count": len(SAMPLES),
        "variant_count": len(VARIANTS),
        "build": "GRCh38",
        "provided_ancestry": "EUR",
        "expected_variant_ids": [f"chr1:{pos}:{ref}:{alt}" for _, pos, ref, alt, _ in VARIANTS],
    }
    (out / "truth.json").write_text(json.dumps(truth, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
