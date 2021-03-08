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
                                args["match"],
                                args["mismatch"],
                                args["gapopen"],
                                args["gapextension"],
                                penalize_end_gaps=args["penalize_end_gaps"],
                                one_alignment_only=True)[0]

    cs = np.cumsum([i!="-" for i in align[0]])

    data = {}
    data["name"] = s["qname"]
    data["seq"] = align[1]
    data["ref"] = align[0]
    data["score"] = align[2]
    target_start = min(np.argwhere(cs==target_start)[0])
    target_end = min(np.argwhere(cs==target_end)[0])
    data["target_ref"] = align[0][target_start:target_end]
    data["target_seq"] = align[1][target_start:target_end]
    data["CB"] = s["CB"]
    data["UB"] = s["UB"]
    return data

def parse_10x_data(args):
    align_file = pysam.AlignmentFile(args["input"], "rb")
    for seq in align_file:
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
    with open(args["reference"],"r") as f:
        ref = f.readline().strip()
        for line in f:
            region,start,end = line.strip().split(" ")
            if region == "target":
                target_start,target_end = int(start),int(end)

    seqs = parse_10x_data(args)
    if int(args["ncores"]) != 1:
        with Pool(processes=int(args["ncores"])) as pool:
            data = [ret for ret in \
                tqdm(pool.imap_unordered(get_alignment_result,[(s, ref, target_start, target_end, args) \
                for s in seqs]))]
    else:
        data = [ret for ret in \
            tqdm(get_alignment_result((s, ref, target_start, target_end, args)) \
            for s in seqs) ]

    result = pd.DataFrame(data=data)
    result.to_csv(args["output"],sep="\t",index=False)
