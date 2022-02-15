#!/usr/bin/env python3

import sys, os, re
import pysam
import subprocess

def main():

    usage = "\n\n\tusage: {} (filename.bam (or - for stdin))  output.bam\n\n".format(sys.argv[0])

    if len(sys.argv) < 3:
        exit(usage)

    input_bam_file = sys.argv[1]
    output_bam_filename = sys.argv[2]

    if input_bam_file == "-":
        input_bam_file = sys.stdin

    samfile = pysam.AlignmentFile(input_bam_file, "rb")

    outfile = pysam.AlignmentFile(output_bam_filename, "wb", template=samfile)


    reads_excluded = 0
    reads_output = 0
    for read in samfile.fetch():
        if read.is_unmapped or read.is_supplementary or (
            read.is_paired and not read.is_proper_pair):
            reads_excluded += 1
            continue

        reads_output += 1
        outfile.write(read)


    print("{}\n\treads_excluded: {}\n\treads_output: {}\n\n".format(sys.argv[0],
                                                                    reads_excluded,
                                                                    reads_output))

    outfile.close()
    subprocess.check_call("samtools index {}".format(output_bam_filename), shell=True)
    

    sys.exit(0)
    

        
if __name__=='__main__':
    main()


