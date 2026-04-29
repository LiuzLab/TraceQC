#!/usr/bin/env python3
"""
Sequence alignment script for 10x Genomics data in TraceQC.

This script processes BAM files from 10x Genomics single-cell experiments,
extracting sequences and barcode information, then performing global pairwise
sequence alignment against a reference sequence.
"""

import argparse
import numpy as np
import pandas as pd
from Bio import pairwise2
from multiprocessing import Pool
from tqdm import tqdm
import pysam


def get_alignment_result(param):
    """
    Perform alignment for a single 10x read.

    Args:
        param: Tuple containing (read_dict, ref_seq, target_start, target_end, args)

    Returns:
        dict or None: Alignment results with barcodes, or None if alignment fails
    """
    s, ref, target_start, target_end, args = param

    # Perform global alignment
    alignments = pairwise2.align.globalms(
        ref, s["query_sequence"],
        args.match, args.mismatch, args.gapopen, args.gapextension,
        penalize_end_gaps=bool(args.penalize_end_gaps),
        one_alignment_only=True
    )

    if not alignments:
        return None

    align = alignments[0]

    # Calculate cumulative sum for position mapping
    cs = np.cumsum([i != "-" for i in align[0]])

    # Prepare result dictionary
    data = {
        "name": s["qname"],
        "seq": align[1],
        "ref": align[0],
        "score": align[2],
        "CB": s["CB"],
        "UB": s["UB"]
    }

    # Extract target region
    try:
        target_start_idx = np.argwhere(cs == target_start)[0][0]
        target_end_idx = np.argwhere(cs == target_end)[0][0]
        data["target_ref"] = align[0][target_start_idx:target_end_idx + 1]
        data["target_seq"] = align[1][target_start_idx:target_end_idx + 1]
    except IndexError:
        data["target_ref"] = ""
        data["target_seq"] = ""

    return data


def parse_10x_data(bam_file):
    """
    Parse 10x Genomics BAM file to extract unmapped reads with barcodes.

    Args:
        bam_file: Path to BAM file

    Yields:
        dict: Read information including sequence and barcodes
    """
    align_file = pysam.AlignmentFile(bam_file, "rb")
    for seq in align_file.fetch(until_eof=True):
        if seq.is_unmapped and seq.query_sequence:
            s = {
                "query_sequence": seq.query_sequence,
                "qname": seq.qname,
                "CB": seq.get_tag("CB") if seq.has_tag("CB") else "/",
                "UB": seq.get_tag("UB") if seq.has_tag("UB") else "/"
            }
            yield s
    align_file.close()


def alignment(args):
    """
    Main alignment function for 10x data.
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

    # Parse sequences from BAM
    sequences = list(parse_10x_data(args.input))

    # Perform alignments
    ncores = int(args.ncores)
    if ncores > 1:
        with Pool(processes=ncores) as pool:
            data = list(tqdm(
                pool.imap_unordered(
                    get_alignment_result,
                    ((s, ref, target_start, target_end, args) for s in sequences)
                ),
                total=len(sequences),
                desc="Aligning 10x sequences"
            ))
    else:
        data = [
            get_alignment_result((s, ref, target_start, target_end, args))
            for s in tqdm(sequences, desc="Aligning 10x sequences")
        ]

    # Filter out failed alignments and save
    data = [d for d in data if d is not None]
    result = pd.DataFrame(data=data)
    result.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="10x Genomics sequence alignment for TraceQC")
    parser.add_argument('--ref', required=True, help='Path to reference file')
    parser.add_argument('--input', required=True, help='Path to BAM input file')
    parser.add_argument('--ncores', type=int, default=1, help='Number of CPU cores to use')
    parser.add_argument('--output', required=True, help='Path to output TSV file')
    parser.add_argument('--match', type=int, default=2, help='Score for matched characters')
    parser.add_argument('--mismatch', type=int, default=-2, help='Penalty for mismatched characters')
    parser.add_argument('--gapopen', type=int, default=-6, help='Penalty for gap opening')
    parser.add_argument('--gapextension', type=float, default=-0.1, help='Penalty for gap extension')
    parser.add_argument('--penalize_end_gaps', type=int, default=1, help='Whether to penalize end gaps (0 or 1)')

    args = parser.parse_args()
    alignment(args)
