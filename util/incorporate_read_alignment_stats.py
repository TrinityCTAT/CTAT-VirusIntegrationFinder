#!/usr/bin/env python

import sys, os, re
import pysam
from collections import defaultdict
import logging
import argparse
#import pandas as pd
import csv
import statistics as st

logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)



def main():

    parser = argparse.ArgumentParser(description="add alignment stats", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("--supp_reads_bam", required=True, type=str, help="supplemental read alignments")
    parser.add_argument("--vif_full_tsv", required=True, type=str, help="VIF full tsv containing evidence read names")
    parser.add_argument("--output", required=True, type=str, help="output tsv file containing read alignment stats")
    parser.add_argument("--detailed", action='store_true', help='show attribute values for each read instead of mean statistic')
    
    args = parser.parse_args()


    bam_filename = args.supp_reads_bam
    vif_full_tsv = args.vif_full_tsv
    outputfilename = args.output

    samfile = pysam.AlignmentFile(bam_filename, "rb")

    ev_reads = set()


    logger.info("-capturing reads of interest from {}".format(vif_full_tsv))

    fh = open(vif_full_tsv, "rt")
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        readnames = row['readnames'].split(",")
        for readname in readnames:
            ev_reads.add(readname)

    fh.close()

    logger.info("-capturing read alignment stats from bam file: {}".format(bam_filename))
    

    read_to_hit_count = defaultdict(int)
    read_to_min_per_id = dict()
    read_to_max_end_clipping = defaultdict(int)
    read_to_min_anchor_len = dict()

    read_counter = 0

    for read in samfile.fetch():
        read_name = read.query_name

        if read_name not in ev_reads:
            continue
        
        read_counter += 1

        if (read_counter % 10000 == 0):
            logger.info("-processed {} alignments".format(read_counter))

        aligned_bases = len(read.get_aligned_pairs(matches_only=True))
        if read_name in read_to_min_anchor_len:
            read_to_min_anchor_len[read_name] = min(aligned_bases, read_to_min_anchor_len[read_name])
        else:
            read_to_min_anchor_len[read_name] = aligned_bases

        
        NH = read.get_tag('NH')
        read_to_hit_count[read_name] = max(read_to_hit_count[read_name], NH)

        mismatch_count = read.get_tag('NM')
        per_id = 100.0 - ( float(mismatch_count) / aligned_bases * 100.0)

        if read_name in read_to_min_per_id:
            read_to_min_per_id[read_name] = min(read_to_min_per_id[read_name], per_id)
        else:
            read_to_min_per_id[read_name] = per_id
        
        cigar = read.cigarstring

        max_clip = 0
        m =  re.search("^(\d+)S", cigar)
        if m:
            max_clip = int(m.group(1))

        m = re.search("(\d+)S$", cigar)
        if m:
            max_clip = max(max_clip, int(m.group(1)))

        read_to_max_end_clipping[read_name] = max(read_to_max_end_clipping[read_name], max_clip)


        


    logger.info("-generating alignment stats report")

    fh = open(vif_full_tsv, "rt")
    reader = csv.DictReader(fh, delimiter="\t")
    fieldnames = list(reader.fieldnames)
    fieldnames.extend(['hits', 'min_per_id', 'max_end_clipping', 'min_anchor_len'])

    ofh = open(outputfilename, 'wt')
    writer = csv.DictWriter(ofh, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    
    logger.info("-writing outputfile: {}".format(outputfilename))

    for row in reader:
        readnames = row['readnames'].split(",")
        hits = list()
        min_per_ids = list()
        max_end_clipping = list()
        min_anchor_lengths = list()
        
        for readname in readnames:

            if readname not in read_to_hit_count:
                raise RuntimeError("Error, missing hit count for read: {}".format(readname))
            
            hits.append(read_to_hit_count[readname])
            min_per_ids.append(read_to_min_per_id[readname])
            max_end_clipping.append(read_to_max_end_clipping[readname])
            min_anchor_lengths.append(read_to_min_anchor_len[readname])

        if args.detailed:
            row['hits'] = ",".join([str(x) for x in hits])
            row['min_per_id'] = ",".join(["{:.1f}".format(x) for x in min_per_ids])
            row['max_end_clipping'] = ",".join([str(x) for x in max_end_clipping])
            row['min_anchor_len'] = ",".join([str(x) for x in min_anchor_lengths])

        else:
            row['hits'] = "{:.3f}".format(st.mean(hits))
            row['min_per_id'] = "{:.1f}".format(st.mean(min_per_ids))
            row['max_end_clipping'] = "{:.3f}".format(st.mean(max_end_clipping))
            row['min_anchor_len'] = "{:.3f}".format(st.mean(min_anchor_lengths))
    

        writer.writerow(row)

    fh.close()
    ofh.close()
    
    logger.info("-done")

    sys.exit(0)

    

if __name__=='__main__':
    main()
