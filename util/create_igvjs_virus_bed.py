#!/usr/bin/env python3
# encoding: utf-8
if __name__ == "__main__":
    import argparse
    import pandas as pd
    import sys

    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", type=str, required=True, help="prefix.virus_read_counts_summary.tsv")
    parser.add_argument("--output", type=str, required=True, help="Output file")
    parser.add_argument("--num_top_viruses", type=int, required=False, default=None, help="num top viruses")
    args = parser.parse_args()
    df = pd.read_csv(args.summary, sep='\t')  # virus	seqlen	mapped	chim_reads
    df = df[df['mapped'] > 0]
    df.sort_values(by='mapped', ascending=False, inplace=True)

    if (args.num_top_viruses is not None and
        args.num_top_viruses >= 1 and
        len(df) > args.num_top_viruses):

        df = df.head(args.num_top_viruses)
        
    
    df['start'] = 0
    df.to_csv(args.output, index=False, header=False, columns=['virus', 'start', 'seqlen'], sep="\t")

    sys.exit(0)
