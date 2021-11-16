#!/usr/bin/env python

import sys, os, re
import pysam
from collections import defaultdict
import logging
import argparse
import pandas as pd
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
    vif_df = pd.read_csv(vif_full_tsv, sep="\t")
    for _, row in vif_df.iterrows():
        readnames = row['readnames'].split(",")
        for readname in readnames:
            ev_reads.add(readname)


    logger.info("-capturing read alignment stats from bam file: {}".format(bam_filename))
    

    read_to_hit_count = defaultdict(int)
    read_to_max_mismatch_count = defaultdict(int)
    read_to_max_end_clipping = defaultdict(int)

    read_counter = 0

    for read in samfile.fetch():
        read_name = read.query_name

        if read_name not in ev_reads:
            continue
        
        read_counter += 1

        if (read_counter % 10000 == 0):
            logger.info("-processed {} alignments".format(read_counter))
        
        NH = read.get_tag('NH')
        read_to_hit_count[read_name] = max(read_to_hit_count[read_name], NH)

        mismatch_count = read.get_tag('NM')
        read_to_max_mismatch_count[read_name] = max(read_to_max_mismatch_count[read_name], mismatch_count)

        alignment_stats = read.get_tag('SA')

        #print("alignment_stats: {}".format(alignment_stats))
        cigar = alignment_stats.split(",")[3]
        
        #print("cigar: {}".format(cigar))

        max_clip = 0
        if (not read.is_paired) or read.is_read1:
            m =  re.search("^(\d+)[SH]", cigar)
            if m:
                max_clip = int(m.group(1))

        if (not read.is_paired) or read.is_read2:
            m = re.search("(\d+)[SH]$", cigar)
            if m:
                max_clip = max(max_clip, int(m.group(1)))

        read_to_max_end_clipping[read_name] = max(read_to_max_end_clipping[read_name], max_clip)
            


    logger.info("-generating alignment stats report")

    vif_df["hits"] = ""
    vif_df["mismatches"] = ""
    vif_df["max_end_clipping"] = ""

    for i, row in vif_df.iterrows():
        readnames = row['readnames'].split(",")
        hits = list()
        mismatches = list()
        max_end_clipping = list()
        
        for readname in readnames:
            hits.append(read_to_hit_count[readname])
            mismatches.append(read_to_max_mismatch_count[readname])
            max_end_clipping.append(read_to_max_end_clipping[readname])


        if args.detailed:
            vif_df.loc[i, 'hits'] = ",".join([str(x) for x in hits])
            vif_df.loc[i, 'mismatches'] = ",".join([str(x) for x in mismatches])
            vif_df.loc[i, 'max_end_clipping'] = ",".join([str(x) for x in max_end_clipping])

        else:
            vif_df.loc[i, 'hits'] = "{:.3f}".format(st.mean(hits))
            vif_df.loc[i, 'mismatches'] = "{:.3f}".format(st.mean(mismatches))
            vif_df.loc[i, 'max_end_clipping'] = "{:.3f}".format(st.mean(max_end_clipping))

    
    vif_df.drop('readnames', axis=1, inplace=True)

    logger.info("-writing outputfile: {}".format(outputfilename))
    vif_df.to_csv(outputfilename, sep="\t", index=False)

    logger.info("-done")

    sys.exit(0)

    

if __name__=='__main__':
    main()
