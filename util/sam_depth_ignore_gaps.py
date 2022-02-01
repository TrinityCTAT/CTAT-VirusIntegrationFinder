#!/usr/bin/env python3

import sys, os, re
import pysam
from collections import defaultdict

def main():

    usage = "\n\n\tusage: {} filename.bam (or - for stdin)\n\n".format(sys.argv[0])

    if len(sys.argv) < 2:
        exit(usage)

    input_bam_file = sys.argv[1]

    if input_bam_file == "-":
        input_bam_file = sys.stdin

    samfile = pysam.AlignmentFile(input_bam_file, "rb")


    seqarray = defaultdict(int)
    prev_chrom = None
    
    for read in samfile.fetch():
        chrom = samfile.get_reference_name(read.reference_id)
        if prev_chrom is None or prev_chrom != chrom:
            dump_coverage(prev_chrom, seqarray)
            # reinit
            prev_chrom = chrom
            seqarray = defaultdict(int)
            
        
        alignment_blocks = read.get_blocks()
        for alignment_block in alignment_blocks:
            for i in range(alignment_block[0], alignment_block[1] + 1):
                seqarray[str(i)] += 1

    if seqarray:
        dump_coverage(prev_chrom, seqarray)

    sys.exit(0)


def dump_coverage(chrom, seqarray):
    positions = sorted(list([int(i) for i in seqarray.keys()]))
    for pos in positions:
        cov = seqarray[ str(pos) ]
        print(f"{chrom}\t{pos}\t{cov}")

    return
        
if __name__=='__main__':
    main()


