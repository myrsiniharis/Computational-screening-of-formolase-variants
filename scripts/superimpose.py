#!/usr/bin/env python3

import csv, os
import argparse
import hashlib
from Bio.PDB import PDBParser, MMCIFParser, Superimposer, PDBIO, Polypeptide
from Bio.Align import PairwiseAligner
from Bio.SeqUtils import seq1


# -----------------------------
# I/O utilities
# -----------------------------

def load_structure(path, name):
    """Load structure from PDB or mmCIF file."""
    if path.endswith(".cif"):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    return parser.get_structure(name, path)


# -----------------------------
# Structure utilities
# -----------------------------

def get_chain(structure, chain_id=None):
    """Return selected chain or first available chain."""
    model = next(structure.get_models())

    if chain_id:
        return model[chain_id]

    return next(model.get_chains())


def extract_sequence(chain):
    """
    Extract amino acid sequence and residue list from a chain.

    Returns:
        seq (str)
        residues (list)
    """
    seq = []
    residues = []

    for res in chain:
        if Polypeptide.is_aa(res, standard=True):
            try:
                seq.append(seq1(res.get_resname()))
                residues.append(res)
            except KeyError:
                continue

    return "".join(seq), residues


# -----------------------------
# Alignment
# -----------------------------

def align_sequences(seq1, seq2):
    """
    Perform global sequence alignment.

    Returns:
        idx1, idx2 : matched residue indices
        gapped1, gapped2 : aligned sequences with gaps
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"

    alignment = aligner.align(seq1, seq2)[0]

    lines = str(alignment).splitlines()
    gapped1 = lines[0].strip()
    gapped2 = lines[2].strip()

    idx1, idx2 = [], []
    pos1 = pos2 = 0

    for a, b in zip(gapped1, gapped2):
        if a != "-" and b != "-":
            idx1.append(pos1)
            idx2.append(pos2)

        if a != "-":
            pos1 += 1
        if b != "-":
            pos2 += 1

    return idx1, idx2, gapped1, gapped2


# -----------------------------
# Atom extraction
# -----------------------------

def extract_atom_pairs(res1, res2, ca_only):
    """
    Extract matching atom pairs for superposition.

    Returns:
        atoms1, atoms2 (same length)
    """
    atoms1, atoms2 = [], []

    for r1, r2 in zip(res1, res2):

        if ca_only:
            if "CA" in r1 and "CA" in r2:
                atoms1.append(r1["CA"])
                atoms2.append(r2["CA"])
            else:
                raise ValueError(f"CA atom missing in one of the residues: {r1.get_id()} vs {r2.get_id()}")
        else:
            atoms1 += list(r1.get_atoms())
            atoms2 += list(r2.get_atoms())

    if len(atoms1) != len(atoms2):
        raise AssertionError(f"Chains have different number of atoms, cannot be superimposed (consider using --ca option)")

    return atoms1, atoms2


# -----------------------------
# Core
# -----------------------------

def superimpose(
    file1,
    file2,
    chain1=None,
    chain2=None,
    ca_only=False,
    align=True,
):
    """
    Superimpose two structures and return statistics as dict.
    """

    s1 = load_structure(file1, "ref")
    s2 = load_structure(file2, "mob")

    c1 = get_chain(s1, chain1)
    c2 = get_chain(s2, chain2)

    seq1, res1 = extract_sequence(c1)
    seq2, res2 = extract_sequence(c2)

    if len(seq1) != len(seq2):
        ca_only = True

    if not align:
        m = min(len(res1), len(res2))
        M = max(len(res1), len(res2))
        matched1, matched2 = res1[:m], res2[:m]
        gapped1 = seq1 + "-" * (M - len(seq1))
        gapped2 = seq2 + "-" * (M - len(seq2))

    else:
        idx1, idx2, gapped1, gapped2 = align_sequences(seq1, seq2)
        matched1 = [res1[i] for i in idx1]
        matched2 = [res2[i] for i in idx2]

    atoms1, atoms2 = extract_atom_pairs(matched1, matched2, ca_only)

    if not atoms1:
        raise RuntimeError("No atoms available for superposition")

    sup = Superimposer()
    sup.set_atoms(atoms1, atoms2)
    sup.apply(s2.get_atoms())

    return {
        "file1": file1,
        "file2": file2,
        "atoms_file1": len(list(s1.get_atoms())),
        "atoms_file2": len(list(s2.get_atoms())),
        "len_seq1": len(seq1),
        "len_seq2": len(seq2),
        "aligned_residues": len(matched1),
        "atoms_aligned": len(atoms1),
        "rmsd": round(sup.rms, 4),
        "superimposed_pdb": "",
        "aligned_seq1": gapped1,
        "aligned_seq2": gapped2,
    }, s2


# -----------------------------
# CLI
# -----------------------------

def main():
    parser = argparse.ArgumentParser(description="Protein structure superposition")

    parser.add_argument("pdb1", help="First PDB file")
    parser.add_argument("pdb2", help="Second PDB file")
    parser.add_argument("--chain1", help="Chain identifier for first structure")
    parser.add_argument("--chain2", help="Chain identifier for second structure")
    parser.add_argument("-o", "--out_csv", help="Output CSV file (appends if file exists, use '-' for stdout)")
    parser.add_argument("-d","--display", action="store_true", help="Display output in human readable format")
    parser.add_argument("--ca", action="store_true", help="Use only CA atoms. Switched on automatically if sequence lengths differ.")
    parser.add_argument("--no-align", action="store_true", help="Skip sequence alignment, cutoff by shorter length")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--out_pdb", help="Write superimposed structure to given PDB file")
    group.add_argument("--pdb_dir", help="Write superimposed structure to PDB file in given directory (hashed filename, for batch processing)")

    args = parser.parse_args()

    if not args.out_csv and not args.display and not args.out_pdb and not args.pdb_dir:
        parser.error("No output specified. Use -o, -d, --out_pdb, or --pdb_dir.")

    result, superimposed_structure = superimpose(
        args.pdb1,
        args.pdb2,
        chain1=args.chain1,
        chain2=args.chain2,
        ca_only=args.ca,
        align=not args.no_align,
    )

    if args.out_pdb or args.pdb_dir:
        if args.pdb_dir:
            os.makedirs(args.pdb_dir, exist_ok=True)
            hash = hashlib.md5(f"{args.pdb1}_{args.pdb2}".encode()).hexdigest()[:8]
            out_path = os.path.join(args.pdb_dir, f"{hash}.pdb")
        else:
            out_path = args.out_pdb
        io = PDBIO()
        io.set_structure(superimposed_structure)
        io.save(out_path)
        result["superimposed_pdb"] = out_path

    if args.out_csv:
        if args.out_csv == '-':
            writer = csv.DictWriter(os.sys.stdout, fieldnames=result.keys())
            writer.writeheader()
            writer.writerow(result)
        else:
            write_header = not os.path.exists(args.out_csv)
            with open(args.out_csv, "a", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=result.keys())
                if write_header:
                    writer.writeheader()
                writer.writerow(result)
    if args.display:
        if args.out_csv == '-':
            raise ValueError("Cannot use display when outputting CSV to stdout")
        for key, value in result.items():
            print(f"{key:20s}: {value}")
        

if __name__ == "__main__":
    main()