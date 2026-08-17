#!/usr/bin/env python3

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Derive oracle visit and linear-trajectory heritability from simulated longitudinal data")
    parser.add_argument("--longitudinal", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--target-h2", type=float, default=None)
    parser.add_argument("--center-time", type=float, default=0.0)
    return parser.parse_args()


def fit_linear(subframe, value_col, center_time):
    x = subframe["time"].to_numpy(dtype=float) - center_time
    y = subframe[value_col].to_numpy(dtype=float)
    design = np.column_stack([np.ones(len(x)), x])
    beta = np.linalg.lstsq(design, y, rcond=None)[0]
    return beta[0], beta[1]


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    dt = pd.read_csv(args.longitudinal, sep="\t")
    required = {"FID", "IID", "visit", "time", "pheno", "genetic_trajectory"}
    missing = required.difference(dt.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    visit_rows = []
    for visit, sub in dt.groupby("visit", sort=True):
        phen = sub["pheno"].to_numpy(dtype=float)
        genetic_obs = sub["genetic_trajectory"].to_numpy(dtype=float)

        if args.target_h2 is not None:
            genetic_component = np.sqrt(args.target_h2) * genetic_obs
        else:
            genetic_component = genetic_obs

        vg = float(np.var(genetic_component, ddof=1))
        vy = float(np.var(phen, ddof=1))
        visit_rows.append(
            {
                "visit": int(visit),
                "time": float(sub["time"].iloc[0]),
                "n": int(len(sub)),
                "var_pheno": vy,
                "var_genetic_component": vg,
                "oracle_h2": vg / vy if vy > 0 else np.nan,
                "genetic_pheno_cor2": float(np.corrcoef(genetic_component, phen)[0, 1] ** 2) if len(sub) > 1 else np.nan,
            }
        )

    pd.DataFrame(visit_rows).to_csv(out_dir / "oracle_visit_heritability.csv", index=False)

    trajectory_rows = []
    for (fid, iid), sub in dt.groupby(["FID", "IID"], sort=False):
        sub = sub.sort_values("time")
        intercept_y, slope_y = fit_linear(sub, "pheno", args.center_time)
        intercept_g, slope_g = fit_linear(sub, "genetic_trajectory", args.center_time)
        trajectory_rows.append(
            {
                "FID": fid,
                "IID": iid,
                "intercept_t0_pheno": intercept_y,
                "annual_slope_pheno": slope_y,
                "intercept_t0_genetic": intercept_g,
                "annual_slope_genetic": slope_g,
            }
        )

    traj = pd.DataFrame(trajectory_rows)
    if args.target_h2 is not None:
        traj["intercept_t0_genetic_component"] = np.sqrt(args.target_h2) * traj["intercept_t0_genetic"]
        traj["annual_slope_genetic_component"] = np.sqrt(args.target_h2) * traj["annual_slope_genetic"]
    else:
        traj["intercept_t0_genetic_component"] = traj["intercept_t0_genetic"]
        traj["annual_slope_genetic_component"] = traj["annual_slope_genetic"]

    def oracle_row(name_pheno, name_gen):
        vy = float(np.var(traj[name_pheno], ddof=1))
        vg = float(np.var(traj[name_gen], ddof=1))
        return {
            "analysis": name_pheno.replace("_pheno", ""),
            "n": int(len(traj)),
            "var_pheno": vy,
            "var_genetic_component": vg,
            "oracle_h2": vg / vy if vy > 0 else np.nan,
            "genetic_pheno_cor2": float(np.corrcoef(traj[name_gen], traj[name_pheno])[0, 1] ** 2),
        }

    linear_rows = [
        oracle_row("intercept_t0_pheno", "intercept_t0_genetic_component"),
        oracle_row("annual_slope_pheno", "annual_slope_genetic_component"),
    ]
    pd.DataFrame(linear_rows).to_csv(out_dir / "oracle_linear_trajectory_heritability.csv", index=False)

    manifest = {
        "longitudinal": str(Path(args.longitudinal).resolve()),
        "out_dir": str(out_dir.resolve()),
        "center_time": args.center_time,
        "target_h2": args.target_h2,
        "files": [
            "oracle_visit_heritability.csv",
            "oracle_linear_trajectory_heritability.csv",
        ],
    }
    with open(out_dir / "manifest.json", "w") as handle:
        json.dump(manifest, handle, indent=2)

    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
