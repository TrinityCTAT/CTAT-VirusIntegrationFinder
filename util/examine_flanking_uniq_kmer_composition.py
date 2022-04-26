#!/usr/bin/env python3

import sys, os, re
import pandas as pd
import argparse

def main():
    
    parser = argparse.ArgumentParser(description="include fraction unique kmer content metric", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("--vif_tsv", type=str, required=True, help="vif tsv input file")
    parser.add_argument("--output", type=str, required=True, help="vif output including the kmer metrics")

    args = parser.parse_args()


    vif_tsv_filename = args.vif_tsv
    output_filename = args.output

    df = pd.read_csv(vif_tsv_filename, sep="\t")

    df = df.apply(examine_unique_kmer_fraction, axis=1)  # fU = fraction unique

    df.to_csv(output_filename, sep="\t", index=False)

    


    sys.exit(0)




def examine_unique_kmer_fraction(row):

    row['flankA_fU'] = "{:.3f}".format(fraction_unique(row['flankA']))
    row['flankB_fU'] = "{:.3f}".format(fraction_unique(row['flankB']))

    return row


def fraction_unique(nuc_seq):

    K=5
    uniq_kmers = set()
    kmer_count = 0
    nuc_seq = nuc_seq.upper()
    for i in range(0, len(nuc_seq)-K):
        kmer = nuc_seq[i:i+K]
        uniq_kmers.add(kmer)
        kmer_count += 1

    frac_uniq = len(uniq_kmers)/kmer_count

    return frac_uniq
    
    
    

if __name__=='__main__':
    main()
