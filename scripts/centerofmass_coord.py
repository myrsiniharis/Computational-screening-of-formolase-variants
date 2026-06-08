#!/usr/bin/env python3
from pymol import cmd
import glob

# folder with pdbs
pdb_files = glob.glob("output/bestVariant_pdbs/FLS/pdb_files/*.pdb")

# residue to analyze
selection_text = "resi 480"   # example ligand
# or: "resi 123"
# or: "chain A and resi 45"

for pdb in pdb_files:

    obj = pdb.replace(".pdb", "")

    cmd.load(pdb, obj)

    # select residue
    sel_name = "target"
    cmd.select(sel_name, f"{obj} and ({selection_text})")

    # center of mass
    com = cmd.centerofmass(sel_name)

    print(f"{pdb}: COM = {com}")

    cmd.delete(obj)

