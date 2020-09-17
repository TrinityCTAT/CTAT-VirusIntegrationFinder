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
    
    arg_parser = argparse.ArgumentParser(description="counts spanning and split evidence reads\n",
                                         formatter_class=argparse.RawTextHelpFormatter
                                        )
    

    arg_parser.add_argument("--patch_db_bam", type=str, required=True, default=None,
                            help="patch genome aligned bam")

    arg_parser.add_argument("--patch_db_gtf", type=str, required=True, default=None,
                            help="patch gtf")

    arg_parser.add_argument("--output_prefix", type=str, required=True,
                            help="output file prefix for bam and tsv")

    arg_parser.add_argument("--min_anchor", type=int, default=25, required=False,
                            help="minimum num aligned bases on both sides of the breakpoint")

    arg_parser.add_argument("--max_end_clip", type=int, default=10, required=False,
                            help="maximum amount of read end clipping")
    
    args = arg_parser.parse_args()

    bam_file = args.patch_db_bam
    gtf_file = args.patch_db_gtf
    output_prefix = args.output_prefix
    min_anchor = args.min_anchor
    max_end_clip = args.max_end_clip

    output_tsv = output_prefix + ".evidence_counts.tsv"
    output_bam_filename = output_prefix + ".evidence.bam"
    
    ofh_tsv = open(output_tsv, 'wt')

    print("\t".join(["contig", "split", "span", "total"]), file=ofh_tsv)  # report header

    

    contig_readnames_want_set = analyze_bam_n_gtf(bam_file, gtf_file, ofh_tsv, min_anchor, max_end_clip)

    ofh_tsv.close()

    # write the bam file
    logger.info("writing read alignment evidence bam")
    samfile = pysam.AlignmentFile(bam_file, "rb")
    outbam = pysam.AlignmentFile(output_bam_filename, 'wb', template=samfile)
    
    for aligned_read in samfile:
        if aligned_read.is_secondary:
            continue
        
        contig = samfile.get_reference_name(aligned_read.reference_id)
        read_name = aligned_read.query_name
        if read_name in contig_readnames_want_set[contig]:
            outbam.write(aligned_read)

    samfile.close()
    outbam.close()
    

    sys.exit(0)

        
def analyze_bam_n_gtf(bam_file, gtf_file, ofh_tsv, min_anchor, max_end_clip):

    samfile = pysam.AlignmentFile(bam_file, "rb")

    
    contig_to_region_pair_dict = parse_region_pairs_from_gtf(gtf_file)
    
    read_to_contig_and_type = defaultdict(dict)
    contig_readnames_to_anchor_lengths = defaultdict(dict)

    logger.info("examining aligned reads")
    
    for aligned_read in samfile:
        contig = samfile.get_reference_name(aligned_read.reference_id)
        read_name = aligned_read.query_name
        
        brkpt = contig_to_region_pair_dict[contig][0]['rend']

        align_start = aligned_read.reference_start
        mate_start = aligned_read.next_reference_start

        align_blocks = aligned_read.get_blocks()
        align_rend = align_blocks[-1][1]

        max_rend = max(mate_start, align_rend)

        if excessive_clipping(aligned_read, max_end_clip):
            continue

        if read_name in read_to_contig_and_type[contig] and read_to_contig_and_type[contig][read_name] == 'split':
            # already handled the upstream mate
            update_anchor_lengths(contig_readnames_to_anchor_lengths[contig][read_name], brkpt, align_blocks)
            continue
        
        # check fragment span
        if read_name in contig_readnames_to_anchor_lengths[contig] or (align_start < brkpt and max_rend > brkpt):
            # fragment overlaps breakpoint.


            if read_name not in contig_readnames_to_anchor_lengths[contig]:
                # init it
                contig_readnames_to_anchor_lengths[contig][read_name] = [0,0]
                read_to_contig_and_type[contig][read_name] = 'span' # default
            
            update_anchor_lengths(contig_readnames_to_anchor_lengths[contig][read_name], brkpt, align_blocks)
            
            # see if we have evidence for a split read
            # which involves a single read spanning the breakpoint
            if align_start < brkpt and align_rend > brkpt:
                # upgrade to split read
                read_to_contig_and_type[contig][read_name] = 'split'
                
    
    logger.info("counting up passing reads")
    # count up the passing reads.
    contig_to_passing_reads = defaultdict(set)
    for contig in contig_readnames_to_anchor_lengths:
        counter = 0
        for readname in contig_readnames_to_anchor_lengths[contig]:
            anchor_lengths = contig_readnames_to_anchor_lengths[contig][readname]
            if anchor_lengths[0] >= min_anchor and anchor_lengths[1] >= min_anchor:
                contig_to_passing_reads[contig].add(readname)
                counter += 1
                #print("{}\t{}\t{}\t{}".format(counter, contig, readname, str(anchor_lengths)))
                
    ## report entries
    logger.info("reporting contig evidence counts")
                
    sorted_contigs = sorted(contig_to_passing_reads.keys(), key=lambda x: len(contig_to_passing_reads[x]), reverse=True)


    # prioritize multimappers according to contig support.
    final_passing_contig_reads = defaultdict(set)
    
    read_seen = set()
    for contig in sorted_contigs:
        type_counter = defaultdict(int)
        for read in contig_to_passing_reads[contig]:
            if read not in read_seen:
                type = read_to_contig_and_type[contig][read]
                type_counter[type] += 1
                read_seen.add(read)
                final_passing_contig_reads[contig].add(read)

        num_split = type_counter.get('split', 0)
        num_span = type_counter.get('span', 0)
        num_total = num_split + num_span
    
        print("\t".join([contig, str(num_split), str(num_span), str(num_total)]), file=ofh_tsv)

    return final_passing_contig_reads





def update_anchor_lengths(anchor_len_list, brkpt, aligned_blocks):

    this_read_anchor_len_list = [0,0]

    for align_block in aligned_blocks:
        (lend, rend) = align_block
        if rend < brkpt:
            this_read_anchor_len_list[0] += rend - lend + 1
        elif lend > brkpt:
            this_read_anchor_len_list[1] += rend - lend + 1
        else:
            this_read_anchor_len_list[0] += brkpt - lend
            this_read_anchor_len_list[1] += rend - brkpt

    if this_read_anchor_len_list[0] > anchor_len_list[0]:
        anchor_len_list[0] = this_read_anchor_len_list[0]

    if this_read_anchor_len_list[1] > anchor_len_list[1]:
        anchor_len_list[1] = this_read_anchor_len_list[1]


def excessive_clipping(aligned_read, max_end_clip):

    cigartuples = aligned_read.cigartuples

    #for cigartuple in cigartuples:
    #    print("tuple: {}".format(str(cigartuple)))
    
    # see https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples
    # cigar codes: S=4, H=5
    
    if cigartuples[0][0] in (4, 5) and cigartuples[0][1] > max_end_clip:
        return True

    if len(cigartuples) > 1 and cigartuples[-1][0] in (4, 5) and cigartuples[-1][1] > max_end_clip:
        return True

    return False



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
