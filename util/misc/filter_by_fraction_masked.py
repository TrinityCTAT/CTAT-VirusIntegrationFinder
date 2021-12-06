#!/usr/bin/env python3

import sys, os, re
import pysam


def main():

    usage = "\n\n\tusage: {} target.fasta max_frac_masked\n\n".format(sys.argv[0])

    if len(sys.argv) < 3:
        exit(usage)

    target_fasta_filename = sys.argv[1]
    max_frac_masked = int(sys.argv[2])

    with pysam.FastxFile(target_fasta_filename) as fh:
        for entry in fh:
            seqname = entry.name
            sequence = entry.sequence
            seqlen = len(sequence)
            num_Ns = count_Ns(sequence)

            frac_masked = "{:.3f}".format(num_Ns / seqlen)

            print("\t".join([seqname, frac_masked]))

    sys.exit(0)

    


def count_Ns(sequence):

    Ncount = 0

    for char in sequence:
        if char == "N" or char == 'n':
            Ncount += 1

    return Ncount



if __name__=='__main__':
    main()

