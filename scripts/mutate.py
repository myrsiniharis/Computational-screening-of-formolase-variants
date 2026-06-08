#!/usr/bin/env python3

import sys, os, re, json, csv, hashlib, argparse
import multiprocessing as mp
from functools import partial
from datetime import datetime
from typing import List, Union

from pyrosetta import init, rosetta, pose_from_file, Pose, create_score_function
from pyrosetta.toolbox import mutate_residue
from pyrosetta.rosetta.protocols.relax import FastRelax
# Developed and tested with pyrosetta 2025.51+release.612b6ef9e9
# written by U
AA = "ACDEFGHIKLMNPQRSTVWY"
mut_pattern = re.compile(rf"(?:[A-Z]\:)?[{AA}]\d+[{AA}]")


def init_worker(structure_file: str, seed: Union[int, None] = None, log_file: Union[str, None] = None):
    """
    Worker initializer for multiprocessing. Initializes PyRosetta and loads the wildtype pose and scoring function.
    Args:
        structure_file: Path to the input PDB/CIF file
        seed: Optional random seed for reproducibility (default: None)
        log_file: Optional file for PyRosetta logs (default: None)
    """
    global wt_pose, scorefxn, relaxor

    init_args = f"-in:file:fullatom -ignore_unrecognized_res"
    if seed: #NOTE: Setting seed makes process deterministic (expected) and consistent (biased)
        # init_args += " -constant_seed"
        init_args += f" -run:seed_offset {seed}"
    if log_file:
        init_args += f" -out:file:log {log_file}_worker{mp.current_process().name}.log"
    else:
        init_args += " -mute all"  # Suppress PyRosetta logs if no log file is specified
    
    with open(os.devnull, "w") as fnull: # Suppress PyRosetta stdout during initialization (banner spaming)
        old_stdout = sys.stdout
        try:
            sys.stdout = fnull
            init(init_args)
        finally:
            sys.stdout = old_stdout

    wt_pose = pose_from_file(structure_file)
    scorefxn = create_score_function("ref2015")
    relaxor = FastRelax(scorefxn)


def listener(queue, output_file: Union[str, None]):
    """
    Listener function for multiprocessing. Listens for results from worker processes and writes them to a CSV file or stdout.
    Args:
        queue: Multiprocessing queue to receive results from worker processes
        output_file: Optional path to output CSV file (default: None, writes to stdout)
    """    
    fieldnames = ["mutation_set", "replicate", "score", "starttime", "endtime", "pdb_file"]

    if output_file:
        f = open(output_file, "w", newline="")
        should_close = True
    else:
        f = sys.stdout
        should_close = False

    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    f.flush()

    try:
        while True:
            item = queue.get()
            if item == "kill":
                break
            writer.writerow(item)
            f.flush()

            # Only fsync real files (not stdout)
            if output_file:
                os.fsync(f.fileno())

    finally:
        if should_close:
            f.close()

def apply_mutations(pose: Pose, mutation_set: List[str]) -> Pose:
    """
    Applies a set of mutations to a given pose and returns the mutated pose.
    Args:
        pose: The original Pose object to mutate
        mutation_set: A list of mutation strings (e.g. ["A23V", "L45M"])
    Returns:
        A new Pose object with the specified mutations applied
    """
    mutant = pose.clone()
    pdb_info = mutant.pdb_info()

    for mutation in mutation_set:
        #TODO: Improve mutation parsing (quite fragile at the moment, not validated againt multi-chain pdb)
        wt_aa = mutation[0]
        new_aa = mutation[-1]

        # Handle optional chain (A:23V)
        if ":" in mutation:
            chain, rest = mutation.split(":")
            resnum = int(rest[1:-1])
            pos = pdb_info.pdb2pose(chain, resnum)
        else:
            pos = int(mutation[1:-1])

        if pos == 0:
            raise ValueError(f"Residue not found in pose: {mutation}")

        actual_aa = mutant.residue(pos).name1()
        if actual_aa != wt_aa:
            raise ValueError(
                f"{mutation}: WT mismatch (pose {pos}) expected {wt_aa}, found {actual_aa}"
            )

        mutate_residue(mutant, pos, new_aa)

    return mutant


def mutate_relax_score(task, queue, out_pdb: str, seed: int = 42):
    """
    Worker function to apply mutations, relax the structure, and score it. Results are sent back to the listener via a queue.
    Args:
        task: A tuple containing (mutation_set, replicate_number)
        queue: Multiprocessing queue to send results back to the listener
        out_pdb: Optional output directory for mutated PDB files
        seed: Random seed for reproducibility (default: 42)
    """
    mutation_set, replicate = task
    if not isinstance(mutation_set, list):
        raise ValueError(f"Invalid mutation set: {mutation_set}. Expected a list of mutation strings.")

    # Set new seed for each replicate to ensure variability (still deterministic due to constant_seed)
    # rosetta.basic.random.set_seed(seed + rep)

    label = "+".join(mutation_set) if mutation_set else "WT"

    try:
        starttime = datetime.now()
        mutant = apply_mutations(wt_pose, mutation_set) if mutation_set else wt_pose.clone()
        relaxor.apply(mutant)
        score = scorefxn(mutant)
        endtime = datetime.now()

        if out_pdb:
            seq = mutant.sequence()
            seq_hash = hashlib.md5(seq.encode()).hexdigest()[:8]
            fname = os.path.join(out_pdb, f"{seq_hash}_rep{replicate}.pdb")
            mutant.dump_pdb(fname)  #TODO: Add remark with sequence and/or mutation set. Maybe score as well?

        queue.put({
            "mutation_set": label,
            "replicate": replicate,
            "score": f"{score:.4f}",
            "starttime": starttime.isoformat(),
            "endtime": endtime.isoformat(),
            "pdb_file": fname if out_pdb else "File not saved"
        })

    except Exception as e:
        queue.put({
            "mutation_set": label,
            "replicate": replicate,
            "score": f"ERROR: {e}",
            "starttime": datetime.now().isoformat(),
            "endtime": datetime.now().isoformat(),
            "pdb_file": "File not saved"
        })



def main():

    parser = argparse.ArgumentParser(description="Mutate input protein structure using PyRosetta + FastRelax and score it with Ref2015.")
    parser.add_argument("input_structure" , help="Input PDB/CIF file")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-m", "--mutations", nargs="+", help="Mutations list. Use comma-separated values for composite mutations (e.g. I87L A23V,B45L)")
    group.add_argument("--mutation_file", type=argparse.FileType("r"), help="JSON file containing mutation sets. Format: [['A23V', 'B45L'], ['C12D']]")
    parser.add_argument("--include_wildtype", action="store_true", help="Score the wildtype structure along with mutated structures")
    parser.add_argument("--out_csv", help="Output CSV file (default: stdout)")
    parser.add_argument("--out_pdb", help="Optional output folder for PDB files of mutated structures")
    parser.add_argument("-r", "--repeats", type=int, default=10, help="Number of repeats per mutation")
    parser.add_argument("-n", "--n_threads", type=int, default=mp.cpu_count(), help="Number of threads")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--pyro_logs", help="PyRosetta log files prefix (default: None, suggested: '.pyrosetta')") #NOTE: not tested
    args = parser.parse_args()

    # Input validation
    if not os.path.isfile(args.input_structure):
        raise FileNotFoundError(f"Input structure file does not exist: {args.input_structure}")
    if args.n_threads <= 0:
        raise ValueError("Number of threads must be > 0")
    if args.repeats <= 0:
        raise ValueError("Repeats must be > 0")

    # Load mutations
    mutations = [[]] if args.include_wildtype else []
    if args.mutation_file:
        with args.mutation_file as f:
            mutations += json.load(f)
    elif args.mutations:
        mutations += [m.split(",") for m in args.mutations]

    # Validate mutation format
    if not mutations:
        raise ValueError("No mutations provided")
    if not all(isinstance(mutation_set, list) and all(isinstance(single_mut, str) and mut_pattern.fullmatch(single_mut) for single_mut in mutation_set) for mutation_set in mutations):
        raise ValueError("Invalid mutation format")
    #TODO: Remove duplicated mutation sets (e.g. A23V,B45L and B45L,A23V)

    # Create output directory for PDBs if needed
    if args.out_pdb:
        os.makedirs(args.out_pdb, exist_ok=True)

    # Set up multiprocessing
    manager = mp.Manager()
    queue = manager.Queue()

    pool = mp.Pool(
        processes=args.n_threads,
        initializer=init_worker,
        initargs=(args.input_structure, args.seed, args.pyro_logs)
    )

    watcher = pool.apply_async(listener, (queue, args.out_csv))

    worker_func = partial(
        mutate_relax_score,
        queue=queue,
        out_pdb=args.out_pdb,
        seed=args.seed
    )

    # Create list of tasks (mutation set + replicate number)
    tasks = [(mutation_set, replicate) for replicate in range(1, args.repeats + 1) for mutation_set in mutations]

    # Map tasks to worker processes
    try:
        pool.map(worker_func, tasks)

    except KeyboardInterrupt:
        print("Interrupted by user")

    finally:
        for _ in range(args.n_threads):
            queue.put("kill")
        pool.close()
        pool.join()

if __name__ == "__main__":
    main()