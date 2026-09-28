#!/usr/bin/env python3
"""Generate small PLINK1 BED/BIM/FAM toy sets for the single-ancestry PRS smoke.

Uses REAL chr22 HapMap3 variants (rsid/BP/alleles) pulled from the PRS-CSx LD
reference snpinfo file so that PRS-CS -- which requires the study SNPs to
overlap the LD reference -- runs end-to-end.

A subset of SNPs is chosen as 'causal' and the FAM phenotype is built from
their genotypes, so the pipeline's real plink2 --glm produces very small
p-values for those SNPs. This guarantees every C+T q-score-range threshold
(0.001, 0.05, ... 0.5 in run_CT.sh) is populated, otherwise run_CT.sh's R step
fails trying to read a missing <range>.profile file.

Also writes an aligned-style sumstats file (SNP CHR BP A1 A2 beta beta_se
P n_eff, matching prepare_sumstats.R output) for direct debugging.

Usage:
    make_smoke_toy.py FILEBASE NSNPS NSAMPLES SEED [IDPREFIX]

FILEBASE e.g. /scratch/.../sim_inputs/AFR/study (writes .bed/.bim/.fam next to
it, plus ss_aligned.txt in the same directory).
IDPREFIX (optional, default "f") is the sample-ID prefix in the FAM; use a
distinct prefix (e.g. "t") so a held-out test sample has unique IIDs.
"""
import os
import random
import sys

SNPINFO = "/projects/standard/gdc/public/prs_methods/ref/ref_PRScsx/1kg_ref/snpinfo_mult_1kg_hm3"


def main():
    filebase = os.path.abspath(sys.argv[1])
    n_snps = int(sys.argv[2])
    n_samples = int(sys.argv[3])
    seed = int(sys.argv[4])
    prefix = sys.argv[5] if len(sys.argv) > 5 else "f"
    rng = random.Random(seed)

    outdir = os.path.dirname(filebase)
    os.makedirs(outdir, exist_ok=True)

    # --- Collect real chr22 SNPs with EUR frequency > 0 from the LD ref ---
    rows = []
    with open(SNPINFO) as fh:
        next(fh)
        for line in fh:
            f = line.strip().split()
            if f[0] == "22" and float(f[8]) > 0:  # FRQ_EUR > 0
                rows.append(f)
            if len(rows) >= n_snps:
                break
    if len(rows) == 0:
        sys.exit("ERROR: could not read any chr22 SNPs from %s" % SNPINFO)
    n_snp = len(rows)
    print("using %d real chr22 HapMap3 SNPs from %s" % (n_snp, SNPINFO))

    # --- genotypes (0/1/2 dosages) and PLINK1 BED (SNP-major magic 6c 1b 01) ---
    genos = [
        [rng.choices([0, 0, 1, 1, 2], k=1)[0] for _ in range(n_samples)]
        for _ in range(n_snp)
    ]
    n_bytes = (n_samples + 3) // 4
    rd = {0: 0b00, 1: 0b10, 2: 0b11}
    buf = bytearray([0x6C, 0x1B, 0x01])
    for row in genos:
        rb = bytearray(n_bytes)
        for i, x in enumerate(row):
            rb[i // 4] |= rd[x] << (2 * (i % 4))
        buf += rb
    with open(filebase + ".bed", "wb") as fh:
        fh.write(bytes(buf))

    # --- BIM (from the real HapMap3 variants) ---
    with open(filebase + ".bim", "w") as fh:
        for f in rows:
            fh.write("22\t%s\t0\t%s\t%s\t%s\n" % (f[1], f[2], f[3], f[4]))

    # --- FAM with a phenotype correlated to a 'causal' SNP subset ---
    # 2 strong causal SNPs (weight ~2.0, small noise) so the glm p-values stay
    # far below 0.001 even with ~50 GWAS samples; everything else is null.
    ncausal = 2
    causal = rng.sample(range(n_snp), ncausal)
    pheno = []
    for i in range(n_samples):
        s = sum(genos[j][i] for j in causal)
        pheno.append(2.0 * s + rng.gauss(0.0, 0.3))
    with open(filebase + ".fam", "w") as fh:
        for i in range(n_samples):
            sex = rng.choice([1, 2])
            fh.write("%s%d\t%s%d\t0\t0\t%d\t%.5f\n" % (prefix, i, prefix, i, sex, pheno[i]))

    # --- Aligned-style sumstats (prepare_sumstats.R schema), for debugging ---
    with open(os.path.join(outdir, "ss_aligned.txt"), "w") as fh:
        fh.write("SNP\tCHR\tBP\tA1\tA2\tbeta\tbeta_se\tP\tn_eff\n")
        for f in rows:
            beta = rng.gauss(0.0, 0.05)
            p = rng.uniform(0.001, 0.5)
            fh.write("%s\t22\t%s\t%s\t%s\t%.6f\t0.05\t%.6g\t100\n" % (f[1], f[2], f[3], f[4], beta, p))

    print("wrote %s.{bed,bim,fam} (%d snps x %d samples, %d causal) + %s/ss_aligned.txt" % (filebase, n_snp, n_samples, ncausal, outdir))


if __name__ == "__main__":
    main()