#!/usr/bin/env python3

import argparse
import random
import subprocess #subprocess is used to run external commands, such as the docking software
import csv
#from modules.docking import vina_box

vina_box = (0, 0, 0, 20, 20, 20)

def get_receptor_bounds(receptor_file, padding=8.0):
    """Extract receptor bounds from PDBQT file"""
    min_x = min_y = min_z = float('inf')
    max_x = max_y = max_z = float('-inf')

    with open(receptor_file) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                parts = line.split()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                min_x, max_x = min(min_x, x), max(max_x, x)
                min_y, max_y = min(min_y, y), max(max_y, y)
                min_z, max_z = min(min_z, z), max(max_z, z)

    return (
    min_x - padding, max_x + padding,
    min_y - padding, max_y + padding,
    min_z - padding, max_z + padding,
)

def random_center_from_bounds(bounds): 
    min_x, max_x, min_y, max_y, min_z, max_z = bounds

    return (
        random.uniform(min_x, max_x),
        random.uniform(min_y, max_y),
        random.uniform(min_z, max_z),
    )

def center_key(center, precision=1.0): #can be changed to 0.5 or 0.1 for finer discretization, but may increase duplicates
    """Discretize center to avoid duplicates"""
    return tuple(int(c / precision) for c in center)

def run_vina(receptor, ligand, center, size, out_file):
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
        # "--log", log_file
    ]

    print("Running:", " ".join(cmd))  # <-- add this line to see command
    subprocess.run(cmd, check=True)


def extract_score(out_file):
    """Extract best (mode 1) score from Vina output"""
    with open(out_file) as f:
        for line in f:
            if line.startswith("REMARK VINA RESULT:"):
                return float(line.split()[3]) # score is 4th element
    return None

def random_docking_batch(receptor, ligand, bounds, n_runs, output_prefix):
    results = []
    used_centers = set()

    _, _, _, sx, sy, sz = vina_box
    size = (sx, sy, sz)

    i = 0
    attempts = 0
    max_attempts = n_runs * 10  # safety to avoid infinite loops

    while i < n_runs and attempts < max_attempts:
        center = random_center_from_bounds(bounds)
        key = center_key(center)

        if key in used_centers:
            attempts += 1
            continue

        used_centers.add(key)

        out_file = f"{output_prefix}_{i}.pdbqt"
        #log_file = f"{output_prefix}_{i}.log"

        try:
            run_vina(receptor, ligand, center, size, out_file)
            score = extract_score(out_file)
        except subprocess.CalledProcessError as e:
            print(f"Vina failed on run {i} at center {center}: {e}")
            score = None

        results.append((i, *center, score))

        print(f"Run {i}: center={center}, score={score}")

        i += 1
        attempts += 1

    return results

def write_csv(results, filename="results.csv"):
    with open(filename, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["run", "center_x", "center_y", "center_z", "score"])
        writer.writerows(results)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Random docking with AutoDock Vina")
    parser.add_argument("--receptor", required=True)
    parser.add_argument("--ligand", required=True)
    parser.add_argument("--runs", type=int, default=20)
    parser.add_argument("--box", nargs=6, type=float, default=vina_box)
    parser.add_argument("--out_prefix", default="random_dock") #prefix for output files, e.g. random_dock_0.pdbqt, random_dock_0.log, etc.
    parser.add_argument("--csv", default="results.csv")
    parser.add_argument("--padding", type=float, default=8.0, help="Padding around protein for random sampling (Å)")

    args = parser.parse_args()

    bounds = get_receptor_bounds(args.receptor, padding=args.padding)

    results = random_docking_batch(
        args.receptor,
        args.ligand,
        bounds,
        args.runs,
        args.out_prefix
    )

    write_csv(results, args.csv)

    print(f"\nResults written to {args.csv}")

    # TODO: Parallelization
    # Idea:
    # - Use multiprocessing.Pool or concurrent.futures.ProcessPoolExecutor
    # - Each worker runs one Vina job
    # - Be careful with file naming to avoid collisions

