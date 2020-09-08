#!/usr/bin/env python3
# encoding: utf-8

import os, re, sys
import argparse
import subprocess
import math
import pysam
from collections import defaultdict
import csv

if sys.version_info[0] < 3:
    raise Exception("This script requires Python 3")


import logging
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s : %(levelname)s : %(message)s',
                    datefmt='%H:%M:%S')



def main():
    
    arg_parser = argparse.ArgumentParser(description="counts spanning and split evidence reads\n" +
                                         "two run modes:\n" +
                                         " - use a single bam and gtf\n" +
                                         " or \n" +
                                         " - give a partitioned data directory and it will iterate through all entries.\n\n",
                                         formatter_class=argparse.RawTextHelpFormatter
                                        )
    

    arg_parser.add_argument("--patch_db_bam", type=str, required=False, default=None,
                            help="patch genome aligned bam")

    arg_parser.add_argument("--patch_db_gtf", type=str, required=False, default=None,
                            help="patch gtf")

    arg_parser.add_argument("--partitioned_dir", type=str, required=False, default=None,
                            help="partitioned data directory, expect to find: 'chim_events_for_eval.tsv' in that dir.")


    args = arg_parser.parse_args()

    bam_file = args.patch_db_bam
    gtf_file = args.patch_db_gtf
    partitioned_dir = args.partitioned_dir

    if not ((bam_file != None and gtf_file != None) ^ (partitioned_dir != None)):
        arg_parser.print_help()
        sys.exit(1)


    
    print("\t".join(["contig", "split", "span", "total"]))  # report header
    
    if partitioned_dir:
        chim_events_file = os.path.join(partitioned_dir, "chim_events_for_eval.tsv")
        if not os.path.exists(chim_events_file):
            raise RuntimeError("Error, cannot locate expected file: {}".format(chim_events_file))
        with open(chim_events_file) as fh:
            csv_reader = csv.DictReader(fh, delimiter="\t")
            for row in csv_reader:
                workdir = row['workdir']
                target_bam = os.path.join(workdir, "Aligned.sortedByCoord.out.bam")
                target_gtf = os.path.join(workdir, "target.gtf")
                analyze_bam_n_gtf(target_bam, target_gtf)
        
    else:
        analyze_bam_n_gtf(bam_file, gtf_file)


    sys.exit(0)

        
def analyze_bam_n_gtf(bam_file, gtf_file):
    
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

    return



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
