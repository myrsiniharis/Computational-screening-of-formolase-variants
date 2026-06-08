#!/usr/bin/env python3

import argparse
import pandas as pd
from Bio import PDB

def parse_csv(csv_file, res_col="residue", val_col="bfactor", use_order=False):
    """Parse the CSV file and return a dictionary mapping residue positions to B-factor values."""
    df = pd.read_csv(csv_file)
    if use_order:
        return list(df[val_col])
    return dict(zip(df[res_col], df[val_col]))

def scale_values(values, scale_range):
    """Scale values to a range of 0-100 based on the provided scale range."""
    min_scale, max_scale = scale_range
    scaled = [(v - min_scale) / (max_scale - min_scale) * 100 for v in values]
    return scaled

def parse_pdb(pdb_file, skip_chain_check):
    """Parse the PDB file and return the structure and residues."""
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    
    chains = [chain for model in structure for chain in model]
    if len(chains) > 1 and not skip_chain_check:
        raise ValueError("Multiple chains detected. This is not implemented yet. Use --skip_chain_check to ignore or ensure chains are identical.")
    
    residues = [residue for chain in chains for residue in chain if residue.get_id()[0] == " "]
    return structure, residues

def apply_bfactors(residues, mapping, use_order, skip_size_check):
    """Apply B-factor values to the PDB structure residues."""
    if use_order:
        if len(mapping) != len(residues):
            raise ValueError("CSV row count does not match residue count in PDB.")
        for residue, value in zip(residues, mapping):
            for atom in residue:
                atom.set_bfactor(value)
    else:
        pdb_residues = {res.get_id()[1] for res in residues}
        csv_residues = set(mapping.keys())
        
        if not skip_size_check and pdb_residues != csv_residues:
            raise ValueError("Residue positions in CSV do not match PDB file.")
        
        for residue in residues:
            res_id = residue.get_id()[1]
            if res_id in mapping:
                for atom in residue:
                    atom.set_bfactor(mapping[res_id])

def write_pdb(structure, output_pdb):
    """Write the modified PDB structure to a file."""
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)

def label_pdb_res(pdb_file, csv_file, output_pdb, res_col="residue", val_col="bfactor", skip_chain_check=False, skip_size_check=False, use_order=False, scale=None):
    """Label PDB file with B-factors from CSV mapping."""
    # Parse CSV file to get B-factor values mapping
    mapping = parse_csv(csv_file, res_col, val_col, use_order)
    
    # Scale values if scaling range is provided
    if scale:
        mapping = scale_values(mapping, scale) if use_order else {k: v for k, v in zip(mapping.keys(), scale_values(mapping.values(), scale))}
    
    # Parse PDB file and extract residues
    structure, residues = parse_pdb(pdb_file, skip_chain_check)
    
    # Apply the new B-factor values
    apply_bfactors(residues, mapping, use_order, skip_size_check)
    
    # Write the modified PDB structure to file
    write_pdb(structure, output_pdb)

def main():
    """Main function to handle argument parsing and coordinate PDB modification steps."""
    parser = argparse.ArgumentParser(description="Update B-factors in a PDB file using a CSV mapping.")
    parser.add_argument("pdb_file", help="Input PDB file")
    parser.add_argument("csv_file", help="CSV file containing residue positions and new B-factor values")
    parser.add_argument("output_pdb", nargs="?", help="Output PDB file with modified B-factors", default="output.pdb")
    parser.add_argument("--res_col", default="residue", help="Column name for residue positions in CSV")
    parser.add_argument("--val_col", default="bfactor", help="Column name for B-factor values in CSV")
    parser.add_argument("--skip_chain_check", action="store_true", help="Skip multi-chain check and use first chain or all identical chains")
    parser.add_argument("--skip_size_check", action="store_true", help="Skip checking if CSV has all residues in PDB")
    parser.add_argument("--use_order", action="store_true", help="Use CSV values in order instead of residue positions (incompatible with --res_col and --skip_size_check)")
    parser.add_argument("--scale", type=float, nargs=2, metavar=('MIN_VAL', 'MAX_VAL'), help="Scale CSV values from given range to 0-100")
    
    args = parser.parse_args()

    label_pdb_res(args.pdb_file, args.csv_file, args.output_pdb, args.res_col, args.val_col, args.skip_chain_check, args.skip_size_check, args.use_order, args.scale)

if __name__ == "__main__":
    main()