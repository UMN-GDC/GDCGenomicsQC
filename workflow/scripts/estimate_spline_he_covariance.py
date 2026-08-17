#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def read_grm(prefix):
    ids = pd.read_csv(f"{prefix}.grm.id", sep=r"\s+", header=None, names=["FID", "IID"], dtype=str)
    n = len(ids)
    packed = np.fromfile(f"{prefix}.grm.bin", dtype=np.float32)
    expected = n * (n + 1) // 2
    if packed.size != expected:
        raise ValueError(f"GRM has {packed.size} values; expected {expected} for {n} samples")
    grm = np.zeros((n, n), dtype=np.float64)
    lower = np.tril_indices(n)
    grm[lower] = packed
    grm[(lower[1], lower[0])] = packed
    return ids, grm


def read_phenotype(path, n_traits):
    cols = ["FID", "IID"] + [f"B{i}" for i in range(1, n_traits + 1)]
    return pd.read_csv(path, sep=r"\s+", header=None, names=cols, dtype={"FID": str, "IID": str})


def read_pc(path, npc):
    cols = ["FID", "IID"] + [f"PC{i}" for i in range(1, npc + 1)]
    return pd.read_csv(path, sep=r"\s+", header=None, usecols=range(npc + 2), names=cols, dtype={"FID": str, "IID": str})


def residualize(y, x):
    beta = np.linalg.lstsq(x, y, rcond=None)[0]
    return y - x @ beta


def he_genetic_covariance(y, grm):
    n, k = y.shape
    count = 0
    sum_x = 0.0
    sum_x2 = 0.0
    sum_z = np.zeros((k, k))
    sum_xz = np.zeros((k, k))

    # Accumulate off-diagonal HE regression moments without storing n(n-1)/2 pairs.
    for i in range(1, n):
        previous = y[:i, :]
        x = grm[i, :i]
        yi = y[i, :]
        sy = previous.sum(axis=0)
        sxy = x @ previous
        sum_z += 0.5 * (np.outer(yi, sy) + np.outer(sy, yi))
        sum_xz += 0.5 * (np.outer(yi, sxy) + np.outer(sxy, yi))
        count += i
        sum_x += x.sum()
        sum_x2 += np.dot(x, x)

    denominator = sum_x2 - (sum_x * sum_x) / count
    if denominator <= 0:
        raise ValueError("HE regression denominator is not positive")
    covariance = (sum_xz - (sum_x / count) * sum_z) / denominator
    return 0.5 * (covariance + covariance.T)


def nearest_psd(matrix):
    values, vectors = np.linalg.eigh(0.5 * (matrix + matrix.T))
    values = np.clip(values, 0, None)
    return (vectors * values) @ vectors.T


def matrix_to_long(matrix, label):
    rows = []
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            rows.append({"matrix": label, "row": f"B{i + 1}", "column": f"B{j + 1}", "value": matrix[i, j]})
    return rows


def main():
    parser = argparse.ArgumentParser(description="Estimate spline genetic covariance with PC-adjusted Haseman-Elston regression")
    parser.add_argument("--grm-prefix", required=True)
    parser.add_argument("--pc", required=True)
    parser.add_argument("--phenotype", required=True)
    parser.add_argument("--basis", required=True)
    parser.add_argument("--basis-df", type=int, required=True)
    parser.add_argument("--npc", type=int, default=10)
    parser.add_argument("--covariance-out", required=True)
    parser.add_argument("--trajectory-out", required=True)
    args = parser.parse_args()

    ids, grm = read_grm(args.grm_prefix)
    phen = read_phenotype(args.phenotype, args.basis_df)
    pcs = read_pc(args.pc, args.npc)

    frame = ids.assign(grm_order=np.arange(len(ids))).merge(phen, on=["FID", "IID"], how="inner")
    frame = frame.merge(pcs, on=["FID", "IID"], how="inner").sort_values("grm_order")
    index = frame["grm_order"].to_numpy(dtype=int)
    grm = grm[np.ix_(index, index)]

    y = frame[[f"B{i}" for i in range(1, args.basis_df + 1)]].to_numpy(dtype=float)
    x = np.column_stack([np.ones(len(frame)), frame[[f"PC{i}" for i in range(1, args.npc + 1)]].to_numpy(dtype=float)])
    y = residualize(y, x)

    sigma_g_raw = he_genetic_covariance(y, grm)
    sigma_total_raw = np.cov(y, rowvar=False, ddof=1)
    sigma_e_raw = sigma_total_raw - np.mean(np.diag(grm)) * sigma_g_raw
    sigma_g = nearest_psd(sigma_g_raw)
    sigma_e = nearest_psd(sigma_e_raw)

    covariance_rows = []
    covariance_rows += matrix_to_long(sigma_g_raw, "genetic_raw")
    covariance_rows += matrix_to_long(sigma_e_raw, "residual_raw")
    covariance_rows += matrix_to_long(sigma_g, "genetic_psd")
    covariance_rows += matrix_to_long(sigma_e, "residual_psd")
    pd.DataFrame(covariance_rows).to_csv(args.covariance_out, index=False)

    basis = pd.read_csv(args.basis, sep="\t")
    basis_cols = [f"B{i}" for i in range(1, args.basis_df + 1)]
    rows = []
    for _, row in basis.iterrows():
        b = row[basis_cols].to_numpy(dtype=float)
        vg = float(b @ sigma_g @ b)
        ve = float(b @ sigma_e @ b)
        rows.append({
            "time": row["time"],
            "genetic_variance": vg,
            "residual_variance": ve,
            "heritability": vg / (vg + ve) if vg + ve > 0 else np.nan,
            "n": len(frame),
        })
    pd.DataFrame(rows).to_csv(args.trajectory_out, index=False)

    print(f"Wrote covariance estimates: {args.covariance_out}")
    print(f"Wrote time-varying heritability: {args.trajectory_out}")


if __name__ == "__main__":
    main()
