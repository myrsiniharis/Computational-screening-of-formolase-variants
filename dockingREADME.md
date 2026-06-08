# AutoDock Vina Docking Workflow

## Overview

This repository provides a reproducible workflow for molecular docking using AutoDock Vina. It supports ligand preparation, receptor preparation, single-site docking, batch docking across protein structures, and visualization of docking results.

The pipeline is designed for programmatic execution and scalable batch processing.

---

## Workflow Summary

Ligand preparation → Receptor preparation → Single docking → Batch docking → Post-processing → Visualization

---

## Requirements

### Software
- Conda
- AutoDock Vina
- MGLTools (AutoDockTools)
- Open Babel
- RDKit
- Meeko
- PyMOL (optional)

### Input formats
- Ligands: `.sdf`, `.pdb`, SMILES
- Receptors: `.pdb`, `.cif`
- Outputs: `.pdbqt`, `.csv`, `.json`

---

## Environment Setup

Create required Conda environments:

```bash
conda env create --file autodockmin_env.yml
conda env create --file meeko_babelmin_env.yml
conda env create --file rdkit.yml
conda env create --file mgltoolsmin_env.yml
```

## Ligand Preparation

Ligands can be prepared from multiple sources.

### From PubChem (SDF)
```
wget "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/XXX/record/SDF?record_type=3d" -O ligand.sdf
```
Convert to PDBQT:

```
conda activate mko_obbl
```
```
mk_prepare_ligand.py -i ligand.sdf -o ligand.pdbqt
```
### From SMILES

```
python smiles_to_pdb.py "SMILES_STRING" ligand.pdb --img ligand.png
```
Convert to PDBQT:
```
conda activate mko_obbl
```
```
obabel ligand.pdb -O ligand.pdbqt --addhydrogens --partialcharge gasteiger
```

## Receptor Preparation
### Convert CIF to PDB (if needed)

```
python cif_to_pdb.py -i input.cif -o output.pdb
```
### Download structure
```
wget https://files.rcsb.org/download/XXXX.pdb
```
### Clean structure
```
grep -v "^ANISOU" input.pdb > receptor_clean.pdb
```
### Split alternate conformations (if needed)
```
conda activate mgltools
```
```
python prepare_pdb_split_alt_confs.py -r receptor_clean.pdb
```
### Generate receptor PDBQT
```
prepare_receptor4.py -r receptor_clean.pdb -o receptor.pdbqt -U waters
```

## Single-Site Docking
### Activate environment
```
conda activate autodock
```
### Docking command
```
vina --receptor receptor.pdbqt \
     --ligand ligand.pdbqt \
     --center_x X --center_y Y --center_z Z \
     --size_x SX --size_y SY --size_z SZ \
     --out output.pdbqt
```

### Output
Docking results contain multiple poses ranked by affinity:
```REMARK VINA RESULT: -X.X```
+ More negative values indicate stronger binding affinity
+ A score of **-5** and above is best

## Batch Docking
### Generate docking coordinates
```
python coordinate_from_pdb.py receptor.pdb --all --output centroids.json
```

### Run batch docking
```
conda activate autodock

python batch_docking.py \
  --receptor receptor.pdbqt \
  --ligands ligand1.pdbqt ligand2.pdbqt ligand3.pdbqt \
  --centroids_json centroids.json \
  --out_prefix outputs/batch_run \
  --csv results.csv \
  --box 20 20 20
  ```

## Post-processing
### Plot docking scores
```
python linegraph_batchdocking.py results.csv output.png
```

### Scale and prepare data
```
python prep_file_4_label_scr.py results.csv results_scaled.csv --scale_scores --score_col ligand1_score ligand2_score
```

### Map scores onto structure
```
python label_pdb_res.py receptor.pdb results_scaled.csv output_labeled.pdb \
  --res_col centroid_idx_plus1 \
  --val_col docking_score_scaled
```
## Vizualization
### PyMOL visualization
```
load output_labeled.pdb

spectrum b, blue_white_red
```
+ B-factor column encodes docking scores

## PyMOL Utilities
### Compute active site center of mass
Example:
```
from pymol import cmd

select site, chain A and resi 28+394+419+480
x, y, z = cmd.centerofmass("site")
print(x, y, z)
```

### Notes
+ Multimeric proteins may have active sites spanning multiple chains
+ Monomer docking may miss interface binding sites
+ Visually inspect docking box placement (To do so load the `docked_liagnd_to_receptor.pdbqt` in PyMOL `run vina_box.py`)

### Troubleshooting
#### Common issues
+ Missing parameters (e.g. Mg warnings): usually safe to ignore
+ No alternate conformations: normal for many structures
+ Vina not found: activate correct environment

### Supplementary
#### Run without activating environment
```
conda run -n env_name tool_name
```

Example:
```
conda run -n mko_obbl mk_prepare_ligand.py -i formaldehyde.sdf -o formaldehyde.pdbqt
```
#### Locate scripts
```
find $CONDA_PREFIX -name script_name.py
```
#### tmux
Use ```tmux``` for jobs that are expected to take a long time

#### Help
```
tool_name --help
```
---
### External resources

[Meeko](https://meeko.readthedocs.io/en/develop/lig_prep_basic.html)

[Vina Documentation](https://autodock-vina.readthedocs.io/_/downloads/en/stable/pdf/)

[Conda Cheatsheet](https://docs.conda.io/projects/conda/en/latest/_downloads/843d9e0198f2a193a3484886fa28163c/conda-cheatsheet.pdf)