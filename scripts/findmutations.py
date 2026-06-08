#!/usr/bin/env python3

import json
import argparse

"""
findmutations2.0.py - Compare FASTA sequences and produce JSON for mutate.py
Added in outputs:
    1. JSON list of lists of mutations (mutate.py-compatible).
    2. Optional sidecar JSON mapping reference/target IDs to each mutant.
- Optionally include wildtype (empty mutation set) = useful in the fastrelax scoring
"""

def read_fasta(filename):
    sequences = {}
    header = None
    seq_parts = []
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    sequences[header] = "".join(seq_parts)
                header = line[1:].split()[0]  # Remove '>' and take the first word
                seq_parts = []
            else:
                seq_parts.append(line)
        if header:
            sequences[header] = "".join(seq_parts)
    return sequences


def get_mutations(seq1, seq2):
    """Return a list of mutations comparing seq1 -> seq2"""
    assert len(seq1) == len(seq2), "Sequences must be of the same length."
    mutations = []
    for i, (a, b) in enumerate(zip(seq1, seq2), start=1):
        if a != b:
            mutations.append(f"{a}{i}{b}")
    return mutations


def main():
    parser = argparse.ArgumentParser(description="Compare FASTA sequences and output mutations as JSON.")
    parser.add_argument("fasta1", help="First FASTA file")
    parser.add_argument("id1", nargs="?", default=None, help="Sequence ID from first FASTA, will be considered the Reference (default: first in FASTA1)")
    parser.add_argument("fasta2", nargs="?", default=None, help="Second FASTA file (default: same as first)")
    parser.add_argument("id2", nargs="?", default=None, help="Target sequence ID or ALL (default: all sequences in same file)")
    parser.add_argument("-o", "--output", help="Output JSON file (list of lists)")
    parser.add_argument("-m", "--metadata", help="Optional JSON mapping reference/target IDs (dictionaries)")
    parser.add_argument("--include_wildtype", action="store_true", help="Include wildtype as empty mutation set")
    args = parser.parse_args()

    # Determine FASTA files
    fasta2_file = args.fasta2 if args.fasta2 else args.fasta1

    # Read sequences
    sequences1 = read_fasta(args.fasta1)
    sequences2 = read_fasta(fasta2_file)

    # Determine reference sequence
    if args.id1:
        reference_id = args.id1
        if reference_id not in sequences1:
            raise ValueError(f"Reference ID '{reference_id}' not found in {args.fasta1}")
        reference_seq = sequences1[reference_id]
    else:
        reference_id = list(sequences1.keys())[0]
        reference_seq = sequences1[reference_id]

    # Determine target sequences to compare against
    comparison_sequences = []

    if args.id2:
        target_id_arg = args.id2
        if target_id_arg.upper() == "ALL" and fasta2_file == args.fasta1:
            # Compare reference to all others in the same file
            for seq_id, seq_str in sequences2.items():
                if seq_id != reference_id:
                    comparison_sequences.append((seq_id, seq_str))
        else:
            if target_id_arg not in sequences2:
                raise ValueError(f"Target ID '{target_arg}' not found in {fasta2_file}")
            comparison_sequences.append((target_id_arg, sequences2[target_id_arg]))
    else: # Compare reference sequence to all other sequences in same file
        if fasta2_file == args.fasta1:
            for seq_id, seq_str in sequences2.items():
                if seq_id != reference_id:
                    comparison_sequences.append((seq_id, seq_str))
        else:
            # Compare sequences from fasta to sequences in second fasta by order
            reference_ids = list(sequences1.keys())
            target_ids = list(sequences2.keys())
            for ref_id, target_id in zip(reference_ids, target_ids): # could be rewritten as "for _, target_id in zip(reference_ids, target_ids):" because the ref_id is already selected separately so here it is ignored and that is why it is grayed out
                comparison_sequences.append((target_id, sequences2[target_id]))

    # Generate results
    mutations_list = []
    metadata_list = []

    # Include wildtype if requested
    if args.include_wildtype:
        mutations_list.append([])
        if args.metadata:
            metadata_list.append({"reference": reference_id, "target": reference_id, "mutations": []})

    for target_sequence_id, target_sequence_string in comparison_sequences:
        # Skip self comparison if same file
        if target_sequence_id == reference_id and fasta2_file == args.fasta1:
            continue
        muts = get_mutations(reference_seq, target_sequence_string)
        if muts:
            mutations_list.append(muts)
            if args.metadata:
                metadata_list.append({"reference": reference_id, "target": target_sequence_id, "mutations": muts})

    # Write mutate.py JSON output
    if not args.output:
        raise ValueError("Output file for mutations must be specified with -o")
    with open(args.output, "w") as f:
        json.dump(mutations_list, f, indent=2)

    # Optional metadata output
    if args.metadata:
        with open(args.metadata, "w") as f:
            json.dump(metadata_list, f, indent=2)

    print(f"Generated {len(mutations_list)} mutation sets for mutate.py in '{args.output}'")
    if args.metadata:
        print(f"Metadata written to '{args.metadata}'")


if __name__ == "__main__":
    main()
