#!/usr/bin/env python3
# encoding: utf-8
if __name__ == "__main__":
    import argparse
    import pandas as pd

    parser = argparse.ArgumentParser()
    parser.add_argument("summary", help="prefix.virus_read_counts_summary.tsv")
    parser.add_argument("output", help="Output file")
    args = parser.parse_args()
    df = pd.read_csv(args.summary, sep='\t')  # virus	seqlen	mapped	chim_reads
    df['start'] = 0
    df.to_csv(args.output, index=False, header=False, columns=['virus', 'start', 'seqlen'])
