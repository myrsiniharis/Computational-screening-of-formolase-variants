#!/usr/bin/env python3
import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser(description="Process CSV and optionally scale docking scores")
    parser.add_argument("input_file", help="Path to input CSV file")
    parser.add_argument("output_file", help="Path to save updated CSV file")

    parser.add_argument("--overwrite", action="store_true",
                        help="Overwrite centroid_idx instead of creating new column")

    parser.add_argument("--score_col", nargs="+", default=["docking_score"],
                    help="One or more columns containing docking scores")

    parser.add_argument("--scale_scores", action="store_true",
                        help="Create a 0–100 scaled docking score column")
    

    args = parser.parse_args()

    df = pd.read_csv(args.input_file)

    # existing behavior
    if args.overwrite:
        df["centroid_idx"] = df["centroid_idx"] + 1
    else:
        df["centroid_idx_plus1"] = df["centroid_idx"] + 1

    # new scoring logic
    if args.scale_scores:
        for col in args.score_col:
            if col not in df.columns:
                raise ValueError(f"Column '{col}' not found in CSV. Available columns: {list(df.columns)}")

            scores = df[col]

            min_score = float(scores.min())
            shifted = scores - min_score
            max_shifted = float(shifted.max())

            if max_shifted == 0:
                df[f"{col}_scaled"] = 100.0
            else:
                df[f"{col}_scaled"] = (shifted / max_shifted) * 100.0

    df.to_csv(args.output_file, index=False)

if __name__ == "__main__":
    main()