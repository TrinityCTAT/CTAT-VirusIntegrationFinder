#!/usr/bin/env python3
# encoding: utf-8

import os, re, sys
import argparse
import subprocess
import pysam
import csv

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

    arg_parser = argparse.ArgumentParser(description="organizes input reads according to candidate integration event site",
                                         formatter_class=argparse.RawTextHelpFormatter )


    arg_parser.add_argument("--left_fq", type=str, required=True, default=None,
                            help="left (or single) fastq file")

    arg_parser.add_argument("--right_fq", type=str, required=False, default="", #intentionally not None
                            help="right fastq file (optional)")

    arg_parser.add_argument("--chim_events", type=str, required=True,
                            help="candidate events tsv containing supporting read names")

    arg_parser.add_argument("--workdir_base", type=str, required=True,
                            help="base directory for work")

    arg_parser.add_argument("--min_total_reads", type=int, default=2,
                            help="min number of total reads required to pursue")
    

    args_parsed = arg_parser.parse_args()

    left_fq_filename = os.path.abspath(args_parsed.left_fq)
    right_fq_filename = os.path.abspath(args_parsed.right_fq) if args_parsed.right_fq else ""
    chim_events_filename = args_parsed.chim_events
    workdir_base_dir = os.path.abspath(args_parsed.workdir_base)
    min_total_reads = args_parsed.min_total_reads


    if os.path.exists(workdir_base_dir):
        logger.error("Error, workdir_base_dir: [{}] currently exists. Please remove it or specify a new workspace directory to create".format(workdir_base_dir))
        sys.exit(1)

    
    read_to_event_dict, event_info_dict = parse_chim_events(chim_events_filename, min_total_reads)

    event_to_workdir_dict = build_workspace(event_info_dict.keys(), workdir_base_dir)

    logger.info("-extracting evidence reads from {}".format(left_fq_filename))
    capture_event_reads(left_fq_filename, 'reads_R1.fastq', read_to_event_dict, event_to_workdir_dict)

    if right_fq_filename:
        logger.info("-extracting evidence reads from {}".format(right_fq_filename))
        capture_event_reads(right_fq_filename, 'reads_R2.fastq', read_to_event_dict, event_to_workdir_dict)
        

    ## log the summary info in the workspace.
    chim_events_for_eval_filename = os.path.join(workdir_base_dir, "chim_events_for_eval.tsv")
    with open(chim_events_for_eval_filename, 'wt') as ofh:
        print("\t".join(["entry", "chrA", "coordA", "orientA", "chrB", "coordB", "orientB", "workdir"]), file=ofh)
        for event_id, event_info_row in event_info_dict.items():
            workdir = event_to_workdir_dict[event_id]
            print("\t".join([event_info_row['entry'],
                             event_info_row['chrA'],
                             event_info_row['coordA'],
                             event_info_row['orientA'],
                             event_info_row['chrB'],
                             event_info_row['coordB'],
                             event_info_row['orientB'],
                             workdir]), file=ofh)
            

    logger.info("-see {} for targeted sites info".format(chim_events_for_eval_filename))
    
    
    sys.exit(0)


def parse_chim_events(chim_events_filename, min_total_reads):

    read_to_event_dict = dict()
    event_info_dict = dict()

    csv.field_size_limit(int(1e6))
    
    with open(chim_events_filename, 'rt') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            
            if int(row['total_reads']) < min_total_reads:
                continue
            
            event_num = row['entry']
            event_info_dict[event_num] = row
            readnames = row['readnames']
            for readname in readnames.split(","):
                if readname in read_to_event_dict:
                    logger.error("-error, readname {} already assigned to different event".format(readname))
                else:
                    read_to_event_dict[readname] = event_num

    return read_to_event_dict, event_info_dict
    
            
    
def build_workspace(event_nums_list, workdir_base_dir):

    event_to_workdir_dict = dict()

    for event_num in event_nums_list:
        event_workdir = os.path.sep.join([workdir_base_dir, "event_{}".format(event_num)])
        os.makedirs(event_workdir)
        event_to_workdir_dict[event_num] = event_workdir
        
    return event_to_workdir_dict


def capture_event_reads(input_fq_filename, output_fq_filename, read_to_event_dict, event_to_workdir_dict):

    reads_want = read_to_event_dict.copy()

    with pysam.FastxFile(input_fq_filename) as fh:
        for entry in fh:
            readname = entry.name
            #print(readname)
            if readname in reads_want:
                event_num = read_to_event_dict[readname]
                workdir = event_to_workdir_dict[event_num]
                reads_filename = os.path.join(workdir, output_fq_filename)
                with open(reads_filename, 'at') as ofh:
                    print("\n".join(["@" + readname,
                                     entry.sequence,
                                     '+',
                                     entry.quality]), file=ofh)
                    

    return
    
    


if __name__=='__main__':
    main()

    
