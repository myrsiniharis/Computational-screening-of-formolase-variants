#!/usr/bin/env python3

import argparse
from Bio import PDB
import os

def convert_cif_to_pdb(input_path, output_path):
    """Load a CIF file and save it as PDB."""
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Input file not found: {input_path}")
    
    parser = PDB.MMCIFParser()
    structure = parser.get_structure('structure', input_path)

    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(output_path)

    print(f"Successfully converted {input_path} to {output_path}")
    print(f"Structure contains {len(structure)} model(s)")

def main():
    parser = argparse.ArgumentParser(
        description="Convert a CIF file to a PDB file using BioPython."
    )
    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Path to the input CIF file."
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Path to the output PDB file."
    )
    args = parser.parse_args()

    convert_cif_to_pdb(args.input, args.output)

if __name__ == "__main__":
    main()