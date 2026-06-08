#!/usr/bin/env python3
import os
import argparse
import json
import subprocess
import csv
import shutil   # <-- added

vina_box = (0, 0, 0, 20, 20, 20)  # default box size (sx, sy, sz)


def run_vina(receptor, ligand, center, size, out_file):
    """Run AutoDock Vina docking for a single center"""
    cx, cy, cz = center
    sx, sy, sz = size

    cmd = [
        "vina",
        "--receptor", receptor,
        "--ligand", ligand,
        "--center_x", str(cx),
        "--center_y", str(cy),
        "--center_z", str(cz),
        "--size_x", str(sx),
        "--size_y", str(sy),
        "--size_z", str(sz),
        "--out", out_file,
    ]

    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def extract_score(out_file):
    """Extract best (mode 1) score from Vina output"""
    with open(out_file) as f:
        for line in f:
            if line.startswith("REMARK VINA RESULT:"):
                return float(line.split()[3])
    return None


def dock_with_repeats_multi(receptor, ligands, centroids, repeats_per_centroid, output_prefix, size):
    """
    Dock multiple ligands on each centroid multiple times.
    Returns list of dicts:
    {
        centroid_idx, repeat_idx, cx, cy, cz,
        ligand1: score, ligand2: score, ...
    }
    """
    results = []

    for c_idx, item in enumerate(centroids):
        center = item["center"]
        resid = item["resid"]
        chain = item["chain"]
        for r_idx in range(repeats_per_centroid):
            row = {
                "centroid_idx": c_idx,
                "repeat_idx": r_idx,
                "resid": resid,
                "chain": chain,
                "center_x": center[0],
                "center_y": center[1],
                "center_z": center[2],
            }

            for ligand in ligands:
                ligand_name = os.path.splitext(os.path.basename(ligand))[0] + "_score"
                out_file = f"{output_prefix}_{ligand_name}_C{c_idx}_R{r_idx}.pdbqt"

                try:
                    run_vina(receptor, ligand, center, size, out_file)
                    score = extract_score(out_file)
                except subprocess.CalledProcessError as e:
                    print(f"Vina failed: ligand={ligand_name}, centroid={c_idx}, repeat={r_idx}")
                    score = None

                row[ligand_name] = score
                print(f"{ligand_name} | Centroid {c_idx}, Repeat {r_idx}: score={score}")

            results.append(row)

    return results


def write_csv(results, filename="results.csv"):
    if not results:
        return
    column_names = list(results[0].keys())
    with open(filename, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=column_names)
        writer.writeheader()
        writer.writerows(results)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Dock ligand at residue centroids with repeats")
    parser.add_argument("--receptor", required=True)
    parser.add_argument("--ligands", nargs="+", required=True, help="List of ligand files")
    parser.add_argument("--centroids_json", required=True, help="JSON file with residue centroid coordinates")
    parser.add_argument("--residues", nargs=2, type=int, default=None, help="Residue range: start end (inclusive), e.g. 1 10")
    parser.add_argument("--repeats_per_centroid", type=int, default=1, help="Number of repeats per centroid")
    parser.add_argument("--out_prefix", default="centroid_dock")
    parser.add_argument("--csv", default="results.csv")
    parser.add_argument("--box", nargs=3, type=float, default=vina_box[3:], help="Docking box size (sx sy sz)")

    args = parser.parse_args()

    # -------------------- ADDED PRE-FLIGHT CHECKS --------------------
    print("Checking inputs...")

    # Check vina exists
    if shutil.which("vina") is None:
        raise FileNotFoundError("AutoDock Vina executable 'vina' not found in PATH.")

    # Check receptor
    if not os.path.isfile(args.receptor):
        raise FileNotFoundError(f"Receptor file not found: {args.receptor}")

    # Check ligands
    for lig in args.ligands:
        if not os.path.isfile(lig):
            raise FileNotFoundError(f"Ligand file not found: {lig}")

    # Check centroids JSON
    if not os.path.isfile(args.centroids_json):
        raise FileNotFoundError(f"Centroids JSON not found: {args.centroids_json}")
    
    # Check output directory for prefix
    out_dir = os.path.dirname(args.out_prefix)
    if out_dir and not os.path.isdir(out_dir):
        raise FileNotFoundError(f"Output directory for prefix does not exist: {out_dir}")

    # Check output directory for CSV
    csv_dir = os.path.dirname(args.csv)
    if csv_dir and not os.path.isdir(csv_dir):
        raise FileNotFoundError(f"CSV output directory does not exist: {csv_dir}")

    print("All input files verified. Starting docking...")
    # ----------------------------------------------------------------

    with open(args.centroids_json) as f:
        raw_centroids = json.load(f)

    centroids = []
    for item in raw_centroids:
        resid = item["resid"]

        if args.residues:
            start, end = args.residues
            if not (start <= resid <= end):
                continue

        centroids.append({
            "resid": resid,
            "chain": item["chain"],
            "center": item["centroid"]
        })

    size = tuple(args.box)
    results = dock_with_repeats_multi(
        args.receptor,
        args.ligands,
        centroids,
        args.repeats_per_centroid,
        args.out_prefix,
        size
    )
    write_csv(results, args.csv)
    print(f"\nResults written to {args.csv}")