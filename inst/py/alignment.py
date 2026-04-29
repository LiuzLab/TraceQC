#!/usr/bin/env python3
"""
Sequence alignment script for TraceQC.

This script performs global pairwise sequence alignment of FASTQ reads against a reference
sequence using Biopython. It supports parallel processing and outputs alignment results
in TSV format for downstream analysis in R.
"""

import os
import argparse
import numpy as np
import pandas as pd
from Bio import pairwise2, SeqIO
from multiprocessing import Pool
from tqdm import tqdm
import gzip


def get_alignment_result(param):
    """
    Perform alignment for a single read and extract target region.

    Args:
        param: Tuple containing (seq_record, ref_seq, target_start, target_end, args)

    Returns:
        dict: Alignment results with target region extracted
    """
    s, ref, target_start, target_end, args = param

    # Handle sequence data (bytes vs str)
    seq = s.seq._data
    if isinstance(seq, bytes):
        seq = seq.decode()

    # Perform global alignment using Biopython
    align = pairwise2.align.globalms(
        ref, seq,
        args.match, args.mismatch, args.gapopen, args.gapextension,
        penalize_end_gaps=bool(args.penalize_end_gaps),
        one_alignment_only=True
    )[0]

    # Calculate cumulative sum to find positions excluding gaps
    cs = np.cumsum([i != "-" for i in align[0]])

    # Prepare result dictionary
    data = {
        "name": s.name,
        "description": s.description,
        "seq": align[1],  # Aligned sequence
        "ref": align[0],  # Aligned reference
        "score": align[2]  # Alignment score
    }

    # Find target region boundaries in aligned sequences
    try:
        target_start_idx = np.argwhere(cs == target_start)[0][0]
        target_end_idx = np.argwhere(cs == target_end)[0][0]
        data["target_ref"] = align[0][target_start_idx:target_end_idx + 1]
        data["target_seq"] = align[1][target_start_idx:target_end_idx + 1]
    except IndexError:
        # If target region not found, use empty strings
        data["target_ref"] = ""
        data["target_seq"] = ""

    return data


def alignment(args):
    """
    Main alignment function.

    Reads reference file, parses FASTQ input, performs alignments,
    and writes results to output file.
    """
    # Read reference sequence and target region
    with open(args.ref, "r") as f:
        ref = f.readline().strip()
        target_start = target_end = None
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3 and parts[0] == "target":
                target_start, target_end = int(parts[1]), int(parts[2])
                break

    if target_start is None or target_end is None:
        raise ValueError("Target region not found in reference file")

    # Open input file (supporting gzipped FASTQ)
    if args.input.lower().endswith('.gz'):
        inp_handle = gzip.open(args.input, "rt")
    else:
        inp_handle = open(args.input, "r")

    try:
        # Parse FASTQ sequences
        sequences = list(SeqIO.parse(inp_handle, "fastq"))

        # Perform alignments (parallel or sequential)
        ncores = int(args.ncores)
        if ncores > 1:
            with Pool(processes=ncores) as pool:
                data = list(tqdm(
                    pool.imap_unordered(
                        get_alignment_result,
                        ((s, ref, target_start, target_end, args) for s in sequences)
                    ),
                    total=len(sequences),
                    desc="Aligning sequences"
                ))
        else:
            data = [
                get_alignment_result((s, ref, target_start, target_end, args))
                for s in tqdm(sequences, desc="Aligning sequences")
            ]

        # Convert to DataFrame and save
        result = pd.DataFrame(data=data)
        result.to_csv(args.output, sep="\t", index=False)

    finally:
        inp_handle.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Sequence alignment for TraceQC")
    parser.add_argument('--ref', required=True, help='Path to reference file')
    parser.add_argument('--input', required=True, help='Path to FASTQ input file')
    parser.add_argument('--ncores', type=int, default=1, help='Number of CPU cores to use')
    parser.add_argument('--output', required=True, help='Path to output TSV file')
    parser.add_argument('--match', type=int, default=2, help='Score for matched characters')
    parser.add_argument('--mismatch', type=int, default=-2, help='Penalty for mismatched characters')
    parser.add_argument('--gapopen', type=int, default=-6, help='Penalty for gap opening')
    parser.add_argument('--gapextension', type=float, default=-0.1, help='Penalty for gap extension')
    parser.add_argument('--penalize_end_gaps', type=int, default=1, help='Whether to penalize end gaps (0 or 1)')

    args = parser.parse_args()
    alignment(args)
