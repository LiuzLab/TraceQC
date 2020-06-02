import os
import argparse
import numpy as np
import pandas as pd
from Bio import pairwise2, SeqIO
from multiprocessing import Process, Pool
from tqdm import tqdm
import gzip


def get_score(param):
    s, ref, regions, args = param
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
    data["seq"] = align[1]
    data["ref"] = align[0]
    data["score"] = align[2]
    
    for region in regions:
        start = min(np.argwhere(cs==regions[region][0])[0])
        end = min(np.argwhere(cs==regions[region][1])[0])
        data["%s_ref"%region] = align[0][start:end]
        data["%s_seq"%region] = align[1][start:end]
        
    return data

def alignment(args):
    regions = {}
    with open(args["reference"],"r") as f:
        ref = f.readline().strip()
        for line in f:
            region,start,end = line.strip().split(" ")
            regions[region] = (int(start),int(end))

    data = {"name":[],"seq":[],"ref":[],"score":[]}
    for region in regions:
        data["%s_seq"%region] = []
        data["%s_ref"%region] = []

    inp_handle = None
    if args["input"].lower().endswith("gz"):
        inp_handle = gzip.open(args["input"], "rt")
    else:
        inp_handle = open(args["input"], "r")

    nreads = 0
    
    if int(args["ncores"]) != 1:
        with Pool(processes=int(args["ncores"])) as pool:
            data = [ret for ret in \
                tqdm(pool.imap_unordered(get_score,((s, ref, regions, args) \
                for s in SeqIO.parse(inp_handle, "fastq")))) ]
    else: 
        data = [ret for ret in \
            tqdm(get_score((s, ref, regions, args)) \
            for s in SeqIO.parse(inp_handle, "fastq")) ]
    
    
    if inp_handle != None:
        inp_handle.close()

    result = pd.DataFrame(data=data)
    result.to_csv(args["output"],sep="\t",index=False)
