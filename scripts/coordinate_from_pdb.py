#!/usr/bin/env python3

import argparse
import json
from collections import defaultdict

amino_acids = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS",
    "ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP",
    "TYR","VAL"
}

def is_amino_acid(resname: str) -> bool:
    return resname.strip() in amino_acids


def get_centroids_from_pdb(pdb_file, selected_residues=None, all_residues=False):
    """
    Returns list of residue-level centroid records:
    [
        {
            "chain": "A",
            "resid": 10,
            "resname": "ALA",
            "centroid": [x, y, z]
        },
        ...
    ]
    """

    residue_coords = defaultdict(list)
    residue_meta = {}

    with open(pdb_file, "r") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            resname = line[17:20].strip()
            if not is_amino_acid(resname):
                continue
            chain_id = line[21].strip()
            resid = int(line[22:26])

            # filtering logic
            if not all_residues:
                if selected_residues is None:
                    continue
                if resid not in selected_residues:
                    continue

            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            key = (chain_id, resid)
            residue_coords[key].append((x, y, z))
            residue_meta[key] = {
                "chain": chain_id,
                "resid": resid,
                "resname": resname
            }

    # build output list
    records = []

    for (chain_id, resid) in sorted(residue_coords.keys(), key=lambda x: (x[0], x[1])):
        atoms = residue_coords[(chain_id, resid)]
        n = len(atoms)

        cx = sum(a[0] for a in atoms) / n
        cy = sum(a[1] for a in atoms) / n
        cz = sum(a[2] for a in atoms) / n

        meta = residue_meta[(chain_id, resid)]

        records.append({
            "chain": meta["chain"],
            "resid": meta["resid"],
            "resname": meta["resname"],
            "centroid": [cx, cy, cz]
        })

    return records


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract amino acid centroids from PDB file (Option B format)"
    )
    parser.add_argument("pdb_file", help="Path to PDB file")
    parser.add_argument("--residues", type=str, help="Comma-separated residue IDs (e.g. 10,25,50)")
    parser.add_argument("--all", action="store_true", help="Include all amino acids")
    parser.add_argument("--output", type=str, default="coordinates.json", help="Output JSON file")
    args = parser.parse_args()

    selected_residues = None
    if args.residues:
        selected_residues = [int(r) for r in args.residues.split(",")]

    records = get_centroids_from_pdb(
        args.pdb_file,
        selected_residues=selected_residues,
        all_residues=args.all
    )

    with open(args.output, "w") as f:
        json.dump(records, f, indent=2)

    print(f"Saved {len(records)} residue centroids to {args.output}")