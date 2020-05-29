import os
import argparse
import numpy as np
import pandas as pd
from Bio import pairwise2, SeqIO
from multiprocessing import Process, Pool

# parser = argparse.ArgumentParser()
# parser.add_argument("--input","-f",help="input sequencing file")
# parser.add_argument("--reference","-r",help="reference sequence file")
# parser.add_argument("--match","-a",help="match score for alignment",default=2)
# parser.add_argument("--mismatch","-i",help="mismatch score for alignment",default=-2)
# parser.add_argument("--gapopen","-g",help="gap opening score for alignment",default=-6)
# parser.add_argument("--gapextension","-e",help="gap extension score for alignment",default=-0.1)
# parser.add_argument("--output","-o",help="output file",default="./output/aligned.txt")
# args = parser.parse_args()

def alignment(args):
    _,args["format"] = os.path.splitext(args["input"])
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

    if args["format"] == ".fastq" or args["format"] == ".fasta":
        sequences = SeqIO.parse(args["input"],args["format"][1:])
    elif args["format"] == ".txt":
        df = pd.read_csv(args["input"])
        sequences = df["sequence"].tolist()

    for s in sequences:
        seq = s
        if args["format"] == ".fastq" or args["format"] == ".fasta":
            seq = s.seq._data
            data["name"].append(s.name)

        align = pairwise2.align.globalms(ref,seq,args["match"],args["mismatch"],\
                                     args["gapopen"],args["gapextension"],one_alignment_only=True)[0]
        data["seq"].append(align[1])
        data["ref"].append(align[0])
        data["score"].append(align[2])

        cs = np.cumsum([i!="-" for i in align[0]])
        for region in regions:
            start = min(np.argwhere(cs==regions[region][0])[0])
            end = min(np.argwhere(cs==regions[region][1])[0])
            data["%s_ref"%region].append(align[0][start:end])
            data["%s_seq"%region].append(align[1][start:end])

    if args["format"] == ".txt":
        for c in df.columns:
            if c != "sequence":
                data[c] = df[c]

    result = pd.DataFrame(data=data)
    result.to_csv(args["output"],sep="\t",index=False)

# alignment(args)
