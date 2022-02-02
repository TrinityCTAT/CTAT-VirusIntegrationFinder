#!/usr/bin/env python

import sys, os, re
import pysam
from collections import defaultdict
import logging
import argparse
import pandas as pd
import subprocess
import math


logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)



def main():

    parser = argparse.ArgumentParser(description="add breakpoint entropy stats", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    

    parser.add_argument("--vif_tsv", required=True, type=str, help="VIF full tsv containing evidence read names")
    parser.add_argument("--ref_genome_fasta", required=True, type=str, help="ref genome fasta file")
    parser.add_argument("--viral_genome_fasta", required=True, type=str, help='viral genome fasta file')
    parser.add_argument("--output", required=True, type=str, help="output tsv file containing additional entropy stats")

    
    args = parser.parse_args()

    vif_tsv_filename = args.vif_tsv
    ref_genome_fasta_filename = args.ref_genome_fasta
    viral_genome_fasta_filename = args.viral_genome_fasta
    output_filename = args.output

    ref_genome_fai_filename = ref_genome_fasta_filename + ".fai"
    if not os.path.exists(ref_genome_fai_filename):
        run_cmd("samtools faidx {}".format(ref_genome_fasta_filename))

    viral_genome_fai_filename = viral_genome_fasta_filename + ".fai"
    if not os.path.exists(viral_genome_fai_filename):
        run_cmd("samtools faidx {}".format(viral_genome_fasta_filename))
               
    virus_accs = set(pd.read_csv(viral_genome_fai_filename, sep="\t", header=None)[0].tolist())
            
    
    vif_df = pd.read_csv(vif_tsv_filename, sep="\t")
    print(vif_df.head())

    flank_len = 30
    def compute_entropy (acc, coord, orient, left_or_right_side):
        if left_or_right_side == "left":
            if orient == '+':
                anchor_left = coord - flank_len + 1
                anchor_right = coord
            elif orient == '-':
                anchor_left = coord
                anchor_right = coord + flank_len - 1

            else:
                raise RuntimeError(f"cannot recognize orient {orient}")

        elif left_or_right_side == "right":
            if orient == '+':
                anchor_left = coord
                anchor_right = coord + flank_len -1
            elif orient == '-':
                anchor_left = coord - flank_len + 1
                anchor_right = coord
            else:
                raise RuntimeError(f"cannot recognize orient {orient}")
        else:
            raise RuntimeError(f"cannot recognize left_or_right_side {left_or_right_side}")


        if acc in virus_accs:
            entropy = compute_entropy_seqrange(acc, viral_genome_fasta_filename, anchor_left, anchor_right)
        else:
            entropy = compute_entropy_seqrange(acc, ref_genome_fasta_filename, anchor_left, anchor_right)

        return entropy
            
        
    vif_df["entropyA"] = vif_df.apply(lambda row: compute_entropy(row['chrA'], row['coordA'], row['orientA'], 'left'), axis=1)
    vif_df["entropyB"] = vif_df.apply(lambda row: compute_entropy(row['chrB'], row['coordB'], row['orientB'], 'right'), axis=1)

    # simplify formatting.
    vif_df["entropyA"]  = vif_df["entropyA"].apply(lambda x: "{:.3f}".format(x) )
    vif_df["entropyB"]  = vif_df["entropyB"].apply(lambda x: "{:.3f}".format(x) )
    
    
    logger.info("-writing outputfile: {}".format(output_filename))
    vif_df.to_csv(output_filename, sep="\t", index=False)

    logger.info("-done")

    sys.exit(0)


def compute_entropy_seqrange(acc, fasta_filename, lend, rend):

    cmd = f"samtools faidx {fasta_filename} {acc}:{lend}-{rend}"

    seq_txt = subprocess.check_output(cmd, shell=True).decode()
    seq_txt = "".join(seq_txt.split("\n")[1:])

    char_counter = defaultdict(int)
    for char in seq_txt:
        char_counter[char] += 1

    num_chars = len(seq_txt)
    entropy = 0.0
    for char, count in char_counter.items():
        p = count / num_chars
        entropy += p * math.log2(1/p)
        

    
    return entropy

    

def run_cmd(cmd):
    logger.info("CMD: {}".format(cmd))
    subprocess.check_call(cmd, shell=True)
    return



if __name__=='__main__':
    main()
