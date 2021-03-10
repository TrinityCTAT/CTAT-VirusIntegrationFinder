#!/usr/bin/env python3
# encoding: utf-8

import os, re, sys
import argparse
import pysam
from collections import defaultdict
import subprocess

if sys.version_info[0] < 3:
    raise Exception("This script requires Python 3")

sys.path.insert(
    0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../PyLib"])
    )

import logging

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
    )

def main():

    arg_parser = argparse.ArgumentParser(
        description="partitions genome-aligned bam into host vs. virus vs. fusions",
        formatter_class=argparse.RawTextHelpFormatter,
        )

    arg_parser.add_argument("--bam", type=str, required=True, help="genome-aligned bam")

    arg_parser.add_argument("--virus_only_gtf", type=str, required=True, help="virus-only gtf")

    arg_parser.add_argument("--outdir", type=str, required=True, help="output directory")
    arg_parser.add_argument("--output_prefix", type=str, required=True, help="output prefix for various bam file partitions")

    args = arg_parser.parse_args()

    input_bam_filename = args.bam
    virus_only_gtf = args.virus_only_gtf
    output_dir = args.outdir
    output_prefix = args.output_prefix
    
    # extract virus genome range within imodel
    virus_chrom, virus_lend, virus_rend = get_virus_region(virus_only_gtf)

    # prep input and output bams
    samfile = pysam.AlignmentFile(input_bam_filename, "rb")

    host_reads_bam_file = os.path.join(output_dir, output_prefix + ".host-only.bam")
    host_reads_bam_writer = pysam.AlignmentFile(host_reads_bam_file, "wb", template=samfile)

    virus_reads_bam_file = os.path.join(output_dir, output_prefix + ".virus-only.bam")
    virus_reads_bam_writer = pysam.AlignmentFile(virus_reads_bam_file, "wb", template=samfile)

    fusion_reads_bam_file = os.path.join(output_dir, output_prefix + ".host-virus.fusions.bam")
    fusion_reads_bam_writer = pysam.AlignmentFile(fusion_reads_bam_file, "wb", template=samfile)

    # write reads to corresponding partitions.
    counter = defaultdict(int)
    for read in samfile:
        in_host, in_virus = examine_read_placement(read, virus_chrom, virus_lend, virus_rend)

        if in_host and in_virus:
            fusion_reads_bam_writer.write(read)
            counter['both'] += 1
        elif in_host:
            host_reads_bam_writer.write(read)
            counter['host-only'] += 1
        elif in_virus:
            virus_reads_bam_writer.write(read)
            counter['virus-only'] += 1
        else:
            raise RuntimeError("read is not in host or in virus... not sure here...")

    logger.info("Counted reads: " + str(counter))

    host_reads_bam_writer.close()
    virus_reads_bam_writer.close()
    fusion_reads_bam_writer.close()

    for file in (host_reads_bam_file, virus_reads_bam_file, fusion_reads_bam_file):
        subprocess.check_call("samtools index " + file, shell=True)

    sys.exit(0)

def examine_read_placement(read, virus_chrom, virus_lend, virus_rend):


    in_host = False
    in_virus = False

    for read_block in read.get_blocks():
        lend, rend = read_block
        
        if lend < virus_rend and rend > virus_lend:
            # have overlap with virus.
            in_virus = True

            # see if overlaps host
            if lend < virus_lend or rend > virus_rend:
                in_host = True

        else:
            # no overlap with virus
            in_host = True


    return in_host, in_virus
        
        

        
    
    
    

def get_virus_region(virus_only_gtf):

    with open(virus_only_gtf) as fh:
        for line in fh:
            line = line.rstrip()
            vals = line.split("\t")
            if len(vals) < 8:
                continue
            if (vals[0] == "imodel"
                and
                vals[1] == "virus"
                and
                vals[2] == "region"):

                return "imodel", int(vals[3]), int(vals[4])

    raise RuntimeError("Error, didn't parse virus region coordinates from gtf: {}".format(virus_only_gtf))


                
    
    
if __name__=='__main__':
    main()


