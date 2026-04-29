#!/usr/bin/env python3
"""
Permutation-based alignment score threshold calculation for TraceQC.

This script generates permuted versions of the reference sequence with varying
levels of mutations and aligns them back to the original to establish score
thresholds for distinguishing true alignments from noise.
"""

import random
import pandas as pd
from Bio import pairwise2


def alignment_score_threshold(args):
    """
    Generate alignment score thresholds through permutation testing.

    Args:
        args: Dictionary containing parameters:
            - reference: Path to reference file
            - n: Number of permutations per percentage
            - permutate_percent: List of permutation percentages
            - read_length: Read length for simulation (if not penalizing end gaps)
            - match, mismatch, gapopen, gapextension: Alignment parameters
            - penalize_end_gaps: Whether to penalize end gaps
            - output_file: Path to output TSV file
    """
    # Base pair substitution options (each base can mutate to the other 3)
    bp_substitutions = {
        "A": ["T", "G", "C"],
        "T": ["A", "G", "C"],
        "G": ["A", "T", "C"],
        "C": ["A", "T", "G"]
    }

    # Read reference sequence
    with open(args["reference"], "r") as f:
        ref_raw = f.readline().strip()

    # Initialize results DataFrame
    results_df = pd.DataFrame()
    n_permutations = int(args["n"])

    for permute_percent in args["permutate_percent"]:
        scores = []
        permuted_seqs = []

        for _ in range(n_permutations):
            # Select reference subsequence if not penalizing end gaps
            if int(args["penalize_end_gaps"]):
                ref = ref_raw
            else:
                start_idx = random.randint(0, len(ref_raw) - int(args["read_length"]))
                ref = ref_raw[start_idx:start_idx + int(args["read_length"])]

            # Generate permuted sequence
            seq_length = len(ref)
            permute_positions = random.sample(
                range(seq_length),
                k=int(seq_length * permute_percent)
            )

            permuted = list(ref)
            for pos in permute_positions:
                original_base = permuted[pos]
                # Substitute with a random different base
                permuted[pos] = random.choice(bp_substitutions[original_base])

            permuted_seq = "".join(permuted)
            permuted_seqs.append(permuted_seq)

            # Align permuted sequence back to reference
            alignment = pairwise2.align.globalms(
                ref, permuted_seq,
                args["match"], args["mismatch"], args["gapopen"], args["gapextension"],
                penalize_end_gaps=bool(int(args["penalize_end_gaps"])),
                one_alignment_only=True
            )[0]

            scores.append(alignment[2])

        # Create DataFrame for this permutation percentage
        temp_df = pd.DataFrame({
            "permutated_seq": permuted_seqs,
            "score": scores,
            "permutate_percent": [permute_percent] * n_permutations
        })

        results_df = pd.concat([results_df, temp_df], ignore_index=True)

        # Save intermediate results
        results_df.to_csv(args["output_file"], index=False, sep="\t")
