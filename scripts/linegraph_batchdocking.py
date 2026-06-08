#!/usr/bin/env python3
import matplotlib 
matplotlib.use("Agg")
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def main():

# Load single file
    df = pd.read_csv("output/batch_docking/Cheng4_batch_docking.csv")

# Convert wide → long format
    long_df = df.melt(
        id_vars=["centroid_idx"],
        value_vars=[
            "formaldehyde_score",
            "benzaldehyde_score",
            "glyceraldehyde_score",
            "TPP_score"
        ],
        var_name="dataset",
        value_name="score"
    )

# Clean dataset names (remove "_score")
    long_df["dataset"] = long_df["dataset"].str.replace("_score", "")

# Average over repeats (same as before)
    mean_scores = (
        long_df.groupby(["centroid_idx", "dataset"])["score"]
        .mean()
        .reset_index()
    )

# Sort
    mean_scores = mean_scores.sort_values(["dataset", "centroid_idx"])

# Z-score normalization (UNCHANGED logic)
    mean_scores["normalized"] = (
        mean_scores.groupby("dataset")["score"]
        .transform(lambda x: (x - x.mean()) / x.std())
    )

# Plot (UNCHANGED)
    plt.figure(figsize=(20, 5))
    sns.lineplot(data=mean_scores, x="centroid_idx", y="normalized", hue="dataset")

    plt.xlabel("Residue (Centroid Index)")
    plt.ylabel("Average Docking Score")
    plt.title("Cheng4 Docking Scores Across Runs")

    plt.xticks(range(0, 551, 50))
    plt.grid(True, linestyle="--", alpha=0.5)

    plt.savefig("output/images/Cheng4docking_comparison.png", dpi=300, bbox_inches="tight")

if __name__ == "__main__":
    main()