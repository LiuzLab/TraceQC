import os
import argparse
import numpy as np
import pandas as pd
from Bio import pairwise2, SeqIO
from multiprocessing import Process, Pool
from tqdm import tqdm
import gzip


def get_alignment_result(param):
    s, ref, target_start, target_end, args = param
    seq = s.seq._data
    align = pairwise2.align.globalms(ref,
                                seq,
                                args["match"],
                                args["mismatch"],
                                args["gapopen"],
                                args["gapextension"],
                                one_alignment_only=True)[0]

    cs = np.cumsum([i!="-" for i in align[0]])

    data = {}
    data["name"] = s.name
    data["description"] = s.description
    data["seq"] = align[1]
    data["ref"] = align[0]
    data["score"] = align[2]
    target_start = min(np.argwhere(cs==target_start)[0])
    target_end = min(np.argwhere(cs==target_end)[0])
    data["target_ref"] = align[0][target_start:target_end]
    data["target_seq"] = align[1][target_start:target_end]
    return data

def alignment(args):
    with open(args["reference"],"r") as f:
        ref = f.readline().strip()
        for line in f:
            region,start,end = line.strip().split(" ")
            if region == "target":
                target_start,target_end = int(start),int(end)

    inp_handle = None
    if args["input"].lower().endswith("gz"):
        inp_handle = gzip.open(args["input"], "rt")
    else:
        inp_handle = open(args["input"], "r")

    nreads = 0

    if int(args["ncores"]) != 1:
        with Pool(processes=int(args["ncores"])) as pool:
            data = [ret for ret in \
                tqdm(pool.imap_unordered(get_alignment_result,((s, ref, target_start, target_end, args) \
                for s in SeqIO.parse(inp_handle, "fastq")))) ]
    else:
        data = [ret for ret in \
            tqdm(get_alignment_result((s, ref, target_start, target_end, args)) \
            for s in SeqIO.parse(inp_handle, "fastq")) ]


    if inp_handle != None:
        inp_handle.close()

    result = pd.DataFrame(data=data)
    result.to_csv(args["output"],sep="\t",index=False)
