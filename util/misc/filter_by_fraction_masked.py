#!/usr/bin/env python3

import sys, os, re
import pysam
import textwrap


def main():

    usage = "\n\n\tusage: {} target.fasta max_frac_masked\n\n".format(sys.argv[0])

    if len(sys.argv) < 3:
        exit(usage)

    target_fasta_filename = sys.argv[1]
    max_frac_masked = float(sys.argv[2])

    with pysam.FastxFile(target_fasta_filename) as fh:
        for entry in fh:
            seqname = entry.name
            sequence = entry.sequence
            seqlen = len(sequence)
            num_Ns = count_Ns(sequence)

            frac_masked = "{:.3f}".format(num_Ns / seqlen)

            print("\t".join([seqname, frac_masked]), file=sys.stderr)

            if float(frac_masked) <= max_frac_masked:
                print(">{}\n{}".format(seqname, "\n".join(textwrap.wrap(sequence, 60)).rstrip()))
            else:
                print("\t** excluding {}".format(seqname), file=sys.stderr)


    sys.exit(0)

    


def count_Ns(sequence):

    Ncount = 0

    for char in sequence:
        if char == "N" or char == 'n':
            Ncount += 1

    return Ncount



if __name__=='__main__':
    main()

