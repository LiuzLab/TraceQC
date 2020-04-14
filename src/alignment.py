import argparse
import numpy as np
import pandas as pd
from Bio import pairwise2, SeqIO
from multiprocessing import Process, Pool

parser = argparse.ArgumentParser()
parser.add_argument("--input","-f",help="input sequencing file")
parser.add_argument("--format","-m",help="input sequencing file format, default: fastq",default="fastq")
parser.add_argument("--reference","-r",help="reference sequence file")
parser.add_argument("--match","-a",help="match score for alignment",default=2)
parser.add_argument("--mismatch","-i",help="mismatch score for alignment",default=-2)
parser.add_argument("--gapopen","-g",help="gap opening score for alignment",default=-6)
parser.add_argument("--gapextension","-e",help="gap extension score for alignment",default=-0.1)
parser.add_argument("--output","-o",help="output file",default="./aligned.txt")
args = parser.parse_args()

ref_df = pd.read_csv(args.reference,sep=" ")
ref = ref_df["seq"][0]
barcode_start = ref_df["barcode_start"][0]
barcode_end = ref_df["barcode_end"][0]
spacer_start = ref_df["spacer_start"][0]
spacer_end = ref_df["spacer_end"][0]

sequences = SeqIO.parse(args.input,args.format)
seq_count = {}
for s in sequences:
    seq = s.seq._data
    if seq not in seq_count:
        align = pairwise2.align.globalms(ref,seq,args.match,args.mismatch,\
                                         args.gapopen,args.gapextension,one_alignment_only=True)[0]
        seq_count[seq] = {"ref_aligned":align[0],"seq_aligned":align[1],"score":align[2],"count":1}
    else:
        seq_count[seq]["count"] += 1
        
res = {}
for seq in seq_count:
    s = seq_count[seq]
    cs = np.cumsum([i!="-" for i in s["ref_aligned"]])
    bs = min(np.argwhere(cs==barcode_start)[0])
    be = min(np.argwhere(cs==barcode_end)[0])
    ss = min(np.argwhere(cs==spacer_start)[0])
    se = min(np.argwhere(cs==spacer_end)[0])
    
    seq_aligned = s["seq_aligned"][bs:be]
    if seq_aligned in res:
        res[seq_aligned]["count"] += s["count"]
    else:
        ref_aligned = s["ref_aligned"][bs:be]
        values = {"count":s["count"],"ref_aligned":ref_aligned,\
                 "spacer":s["seq_aligned"][ss:se],"score":s["score"]}
        res[seq_aligned] = values
    
keys = list(res.keys())
data = {"seq_aligned":keys,"ref_aligned":[res[k]["ref_aligned"] for k in keys],\
       "count":[res[k]["count"] for k in keys],"spacer":[res[k]["spacer"] for k in keys],\
       "score":[res[k]["score"] for k in keys]}
df = pd.DataFrame(data=data)
df.to_csv(args.output,sep="\t",index=False)
