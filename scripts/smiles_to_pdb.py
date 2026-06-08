#!/usr/bin/env python3

import argparse
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

# TODO: Add SDF batch processing using Chem.SDMolSupplier

def smiles_to_mol(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")
    return mol

def prepare_molecule(mol, add_h=True, optimize=True, charges=True):
    if add_h:
        mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG()) # Generate 3D coordinates
    if optimize:
        AllChem.MMFFOptimizeMolecule(mol)
    if charges:
        AllChem.ComputeGasteigerCharges(mol)
    return mol

def write_pdb(mol, output_file):
    with open(output_file, "w") as f:
        f.write(Chem.MolToPDBBlock(mol))

def generate_image(mol, image_file="molecule.png"):
    """Generate a 2D image of the molecule for verification."""
    img = Draw.MolToImage(mol, size=(300, 300))
    img.save(image_file)
    print(f"2D image saved to {image_file}")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Convert SMILES to docking-ready PDB")
    parser.add_argument("smiles", help="Input SMILES string")
    parser.add_argument("output", help="Output PDB file")
    parser.add_argument("--no-h", action="store_true", help="Do not add hydrogens")
    parser.add_argument("--no-opt", action="store_true", help="Skip geometry optimization")
    parser.add_argument("--no-charge", action="store_true", help="Skip Gasteiger charges")
    parser.add_argument("--img", default="molecule.png", help="Optional 2D verification image file")
    args = parser.parse_args()

    mol = smiles_to_mol(args.smiles)
    mol = prepare_molecule(
        mol,
        add_h=not args.no_h,
        optimize=not args.no_opt,
        charges=not args.no_charge
    )

    write_pdb(mol, args.output)
    generate_image(mol, args.img)

    print(f"Converted SMILES '{args.smiles}' to PDB file '{args.output}'")