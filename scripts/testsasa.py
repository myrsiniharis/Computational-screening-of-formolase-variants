#!/usr/bin/env python3
import os
import pandas as pd
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.Data.IUPACData import protein_letters_3to1

'''
need to check this script might be sketchy 
'''
parser, sr = PDBParser(), ShrakeRupley()

def pdb_to_sasa_df(path_pdb):
    model = parser.get_structure(os.path.splitext(path_pdb)[0], path_pdb)[0]
    sr.compute(model, level="R")
    return pd.DataFrame(
        [(c.id, i,
          "M" if r.get_resname()=="MSE" else protein_letters_3to1.get(r.get_resname(),"X"),
          r.sasa)
         for c in model for i,r in enumerate(c.get_residues())],
        columns=['chain','pos','res','sasa']
    )

pdb_path = "input/AF-P20906-F1-model_v6.pdb"
df = pdb_to_sasa_df(pdb_path)
df.to_csv("sasa_output.csv", index=False)
print(df.head())