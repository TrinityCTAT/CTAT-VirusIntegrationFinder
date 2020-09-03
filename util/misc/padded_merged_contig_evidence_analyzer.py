#!/usr/bin/env python3
# encoding: utf-8

import os, re, sys
import argparse
import subprocess
import math
import pysam
from collections import defaultdict

if sys.version_info[0] < 3:
    raise Exception("This script requires Python 3")

sys.path.insert(0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../PyLib"]))
from Pipeliner import Pipeliner, Command

import logging
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s : %(levelname)s : %(message)s',
                    datefmt='%H:%M:%S')



def main():
    
    arg_parser = argparse.ArgumentParser(description="runs STAR, gathering chimeric junctions",
                                        formatter_class=argparse.RawTextHelpFormatter
                                        )


    arg_parser.add_argument("--patch_db_bam", type=str, required=True,
                            help="patch genome aligned bam")

    arg_parser.add_argument("--patch_db_gtf", type=str, required=True,
                            help="patch gtf")


    args = arg_parser.parse_args()

    bam_file = args.patch_db_bam
    gtf_file = args.patch_db_gtf

    contig_to_region_pair_dict = parse_region_pairs_from_gtf(gtf_file)
    
    samfile = pysam.AlignmentFile(bam_file, "rb")

    read_to_contig_and_type = defaultdict(dict)
    init_contig_support = defaultdict(int)
    
    for aligned_read in samfile:
        contig = samfile.get_reference_name(aligned_read.reference_id)
        read_name = aligned_read.query_name
        
        brkpt = contig_to_region_pair_dict[contig][0]['rend']

        align_start = aligned_read.reference_start
        mate_start = aligned_read.next_reference_start

        align_blocks = aligned_read.get_blocks()
        align_rend = align_blocks[-1][1]

        max_rend = max(mate_start, align_rend)

        # check fragment span
        if read_name in read_to_contig_and_type[contig] and  read_to_contig_and_type[contig][read_name] == 'split':
            # already handled this one
            init_contig_support[contig] +=1
            continue
        
        if align_start < brkpt and max_rend > brkpt:
            # fragment overlaps breakpoint.
            read_to_contig_and_type[contig][read_name] = 'span'
            init_contig_support[contig] +=1
            
            # see if we have evidence for a split read
            # which involves a single read spanning the breakpoint
            if align_rend > brkpt:
                # upgrade to split read
                read_to_contig_and_type[contig][read_name] = 'split'
                

    ## count'em up
    sorted_contigs = sorted(init_contig_support.keys(), key=lambda x: init_contig_support[x], reverse=True)

    print("\t".join(["contig", "split", "span", "total"]))
    
    read_seen = set()
    for contig in sorted_contigs:
        type_counter = defaultdict(int)
        for read in read_to_contig_and_type[contig]:
            if read not in read_seen:
                type = read_to_contig_and_type[contig][read]
                type_counter[type] += 1
                read_seen.add(read)

        num_split = type_counter.get('split', 0)
        num_span = type_counter.get('span', 0)
        num_total = num_split + num_span
    
        print("\t".join([contig, str(num_split), str(num_span), str(num_total)]))
        
            
    sys.exit(0)



def parse_region_pairs_from_gtf(gtf_file):

    contig_to_region_pair_dict = defaultdict(list)
    
    with open(gtf_file, 'rt') as fh:
        for line in fh:
            line = line.rstrip()
            vals = line.split("\t")

            contig = vals[0]
            
            region_dict = {
                'contig' : contig,
                'lend' : int(vals[3]),
                'rend' : int(vals[4]),
                'orient' : vals[6],
                'region_info' : vals[8] }
            
            contig_to_region_pair_dict[contig].append(region_dict)

    return contig_to_region_pair_dict


if __name__=='__main__':
    main()
