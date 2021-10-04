import os
import argparse
import numpy as np
import pandas as pd
from Bio import pairwise2, SeqIO
from multiprocessing import Process, Pool
from tqdm import tqdm
import gzip
import pysam

def get_alignment_result(param):
    s, ref, target_start, target_end, args = param
    align = pairwise2.align.globalms(ref,
                                s["query_sequence"],
                                args.match,
                                args.mismatch,
                                args.gapopen,
                                args.gapextension,
                                penalize_end_gaps=args.penalize_end_gaps,
                                one_alignment_only=True)

    if align:
        align = align[0]
        cs = np.cumsum([i!="-" for i in align[0]])

        data = {}
        data["name"] = s["qname"]
        data["seq"] = align[1]
        data["ref"] = align[0]
        data["score"] = align[2]
        target_start = min(np.argwhere(cs==target_start)[0])
        target_end = min(np.argwhere(cs==target_end)[0])
        data["target_ref"] = align[0][target_start:target_end+1]
        data["target_seq"] = align[1][target_start:target_end+1]
        data["CB"] = s["CB"]
        data["UB"] = s["UB"]
        return data
    return

def parse_10x_data(args):
    align_file = pysam.AlignmentFile(args.input, "rb")
    for seq in align_file.fetch(until_eof=True):
        if seq.is_unmapped: 
            s = {}
            s["query_sequence"] = seq.query_sequence
            s["qname"] = seq.qname
            try:
                s["CB"] = seq.get_tag("CB")
            except:
                s["CB"] = "/"
            try:
                s["UB"] = seq.get_tag("UB")
            except:
                s["UB"] = "/"
            yield s
    
def alignment(args):
    with open(args.ref,"r") as f:
        ref = f.readline().strip()
        for line in f:
            region,start,end = line.strip().split(" ")
            if region == "target":
                target_start,target_end = int(start),int(end)

    seqs = parse_10x_data(args)
    if int(args.ncores) != 1:
        with Pool(processes=int(args.ncores)) as pool:
            data = [ret for ret in \
                tqdm(pool.imap_unordered(get_alignment_result,[(s, ref, target_start, target_end, args) \
                for s in seqs]))]
    else:
        data = [ret for ret in \
            tqdm(get_alignment_result((s, ref, target_start, target_end, args)) \
            for s in seqs)]
    data = [d for d in data if d]
    result = pd.DataFrame(data=data)
    result.to_csv(args.output,sep="\t",index=False)
    
if __name__ == "__main__"   :
    parser = argparse.ArgumentParser()
    parser.add_argument('--ref', help='the path of reference file')
    parser.add_argument('--input', help='the path of bam file')
    parser.add_argument('--ncores', help='the number of cores to use',default=1)
    parser.add_argument('--output', help='the path of output file')
    parser.add_argument('--match', help='the penalty score for matched character', type=int, default=2)
    parser.add_argument('--mismatch', help='the penalty score for mismatched character', type=int, default=-2)
    parser.add_argument('--gapopen', help='the penalty score for gap opening', type=int, default=-6)
    parser.add_argument('--gapextension', help='the penalty score for gap extension', type=float, default=-0.1)
    parser.add_argument('--penalize_end_gaps', help='if penalize the end gaps: 0 or 1', type=int, default=1)
    args = parser.parse_args()
    alignment(args)
