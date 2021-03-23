import random
import pandas as pd
from Bio import pairwise2

def alignment_score_threshold(args):
    bp = {"A":["T","G","C"],
      "T":["A","G","C"],
      "G":["A","T","C"],
      "C":["A","T","G"]}
    
    with open(args["reference"],"r") as f:
        ref_raw = f.readline().strip()
    
    df = pd.DataFrame()
    n = int(args["n"])
    for per in args["permutate_percent"]:
        scores = []
        perms = []
        for _ in range(n):
            if args["penalize_end_gaps"]:
                idx = random.randint(0,len(ref_raw)-args["read_length"])
                ref = ref_raw[idx:idx+args["read_length"]]
            else:
                ref = ref_raw
            length = len(ref)
            perm_id = random.sample(range(length),k=int(length * per))
            perm = [i for i in ref]
            for i in perm_id:
                perm[i] = random.choice(bp[perm[i]])
            perm = "".join(perm)
            perms.append(perm)
            align = pairwise2.align.globalms(ref,
                                        perm,
                                        args["match"],
                                        args["mismatch"],
                                        args["gapopen"],
                                        args["gapextension"],
                                        penalize_end_gaps=args["penalize_end_gaps"],
                                        one_alignment_only=True)[0]
            scores.append(align[2])
        tmp = pd.DataFrame({"permutated_seq":perms,"score":scores,\
            "permutate_percent":[per for _ in range(n)]})
        df = pd.concat([df,tmp])
        df.to_csv(args["output_file"],index=False,sep="\t")
