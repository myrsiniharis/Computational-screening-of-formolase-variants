#!/usr/bin/env python3

import pandas as pd
import argparse
import glob
import os

def sort_csv(file_path, column, output_dir, descending):
    df = pd.read_csv(file_path)

    if column not in df.columns:
        print(f" Column '{column}' not found in {file_path}")
        return

    # Force numeric sorting when possible
    df[column] = pd.to_numeric(df[column], errors="coerce")

    df_sorted = df.sort_values(by=column, ascending=not descending)

    filename = os.path.basename(file_path)
    output_path = os.path.join(output_dir, f"sorted_{filename}")

    df_sorted.to_csv(output_path, index=False)

    print(f" Sorted {filename} → {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Sort CSV files by a column")

    parser.add_argument("--i", required=True,
                        help="Path to a CSV file OR folder containing CSV files")

    parser.add_argument("--column", required=True,
                        help="Column name to sort by")

    parser.add_argument("--o", default="sorted_output",
                        help="Output folder (default: sorted_output)")

    parser.add_argument("--desc", action="store_true",
                        help="Sort in descending order")

    args = parser.parse_args()

    os.makedirs(args.o, exist_ok=True)

    # If folder → process all CSVs
    if os.path.isdir(args.i):
        files = glob.glob(os.path.join(args.i, "*.csv"))
    else:
        files = [args.i]

    if not files:
        print(" No CSV files found.")
        return

    for file_path in files:
        sort_csv(file_path, args.column, args.o, args.desc)


if __name__ == "__main__":
    main()