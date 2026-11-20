import os
import sys
import argparse
import pandas as pd  # type: ignore
import numpy as np   # type: ignore
import subprocess

def main():
    parser = argparse.ArgumentParser(
        description="Split PLINK .bed files into train/val/test subsets"
    )
    parser.add_argument("plink_prefix", help="Path to PLINK file prefix (e.g., /path/to/data/prefix)")
    parser.add_argument("--train", type=float, default=50.0, help="Percent of data for training (default: 50)")
    parser.add_argument("--val", type=float, default=20.0, help="Percent of data for validation (default: 20)")
    parser.add_argument("--test", type=float, default=30.0, help="Percent of data for testing (default: 30)")
    parser.add_argument("--no_plink", action="store_true", help="Only split samples, don't call PLINK")
    parser.add_argument("--seed", type=int, default=42, help="Random seed (default: 42)")
    args = parser.parse_args()

    # Validate split proportions
    total_pct = args.train + args.val + args.test
    if not np.isclose(total_pct, 100.0):
        print(f"ERROR: Percentages must sum to 100. Got {total_pct:.2f}")
        sys.exit(1)

    plink_dir = os.path.dirname(args.plink_prefix)
    if plink_dir == "":
        plink_dir = "."
    prefix_name = os.path.basename(args.plink_prefix)

    fam_path = args.plink_prefix + ".fam"
    if not os.path.exists(fam_path):
        print(f"ERROR: .fam file not found: {fam_path}")
        sys.exit(1)

    # Load and shuffle .fam file
    df = pd.read_csv(fam_path, sep=r'\s+', header=None)
    df = df[[0, 1]].sample(frac=1.0, random_state=args.seed).reset_index(drop=True)

    # Compute sample sizes
    n = len(df)
    n_train = int(n * args.train / 100)
    n_val = int(n * args.val / 100)
    n_test = n - n_train - n_val

    # Split the dataset
    df_train = df.iloc[:n_train]
    df_val = df.iloc[n_train:n_train + n_val]
    df_test = df.iloc[n_train + n_val:]

    # Save split sample lists
    def save_split(df_split, name):
        out_path = os.path.join(plink_dir, f"{name}_samples.txt")
        df_split.to_csv(out_path, sep=' ', index=False, header=False)
        print(f"Saved: {out_path}")
        return out_path

    train_file = save_split(df_train, "train")
    val_file = save_split(df_val, "val")
    test_file = save_split(df_test, "test")

    print(f"\nTotal samples: {n}")
    print(f"Train: {n_train} | Val: {n_val} | Test: {n_test}\n")

    # Optionally run PLINK to generate subsets
    if not args.no_plink:
        def run_plink(input_file, split_name):
            output_prefix = os.path.join(plink_dir, f"{prefix_name}_{split_name}")
            cmd = [
                "plink", "--bfile", args.plink_prefix,
                "--keep", input_file,
                "--make-bed", "--out", output_prefix
            ]
            print(f"Running PLINK for {split_name} subset...")
            subprocess.run(cmd, check=True)

        if n_train > 0:
            run_plink(train_file, "train")
        if n_val > 0:
            run_plink(val_file, "val")
        if n_test > 0:
            run_plink(test_file, "test")

        print("\nPLINK subsets generated successfully.")

if __name__ == "__main__":
    main()
