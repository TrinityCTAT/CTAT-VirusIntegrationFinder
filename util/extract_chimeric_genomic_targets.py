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

    arg_parser = argparse.ArgumentParser(description="extracts target regions for candidate integration event site for draft evaluation",
                                         formatter_class=argparse.RawTextHelpFormatter )


    arg_parser.add_argument("--genome_lib_dir", type=str, default=os.environ.get('CTAT_GENOME_LIB'),
                            help="genome lib directory - see http://FusionFilter.github.io for details.  Uses env var CTAT_GENOME_LIB as default")


    arg_parser.add_argument("--patch_db_fasta", type=str, required=True,
                            help="patch genome database")
    
    arg_parser.add_argument("--output_prefix", type=str, required=True,
                            help="output prefix for fasta and gtf files")

    arg_parser.add_argument("--chim_events", type=str, required=True,
                            help="chimeric events file (ie. prelim.vif.abridged.tsv)")

    arg_parser.add_argument("--pad_region_length", type=int, default=1000,
                            help="length around breakpoint to extract genome sequence")
    
    args_parsed = arg_parser.parse_args()


    genome_lib_dir = args_parsed.genome_lib_dir
    patch_db_fasta = args_parsed.patch_db_fasta
    output_prefix = args_parsed.output_prefix
    pad_region_length = args_parsed.pad_region_length
    chim_events_filename = args_parsed.chim_events

    if not genome_lib_dir:
        logger.error("Error, --genome_lib_dir must be specified");
        sys.exit(1)

    ref_genome_fasta = os.path.join(genome_lib_dir, "ref_genome.fa")

    
    event_info_dict = parse_chim_events(chim_events_filename)

    write_genome_target_regions(event_info_dict, ref_genome_fasta, patch_db_fasta, output_prefix, pad_region_length)
    
    
    sys.exit(0)



def write_genome_target_regions(event_info_dict, ref_genome_fasta, patch_db_fasta, output_prefix, pad_region_length):

    
    patch_db_entries = set()
    with open(patch_db_fasta, 'rt') as fh:
        for line in fh:
            m = re.search("^>(\S+)", line)
            if m:
                acc = m.group(1)
                patch_db_entries.add(acc)


    out_fasta_filename = output_prefix + ".fasta"
    out_gtf_filename = output_prefix + ".gtf"
    

    ofh_fasta = open(out_fasta_filename, 'wt')
    ofh_gtf = open(out_gtf_filename, 'wt')

    for event in event_info_dict.values():

        event_num = event['entry']

        event_acc = "candidate_{}".format(event_num)

        chrA = event['chrA']
        coordA = int(event['coordA'])
        orientA = event['orientA']

        chrB = event['chrB']
        coordB = int(event['coordB'])
        orientB = event['orientB']

        brkpt_type = event['primary_brkpt_type']
        
        chrA_fasta_file, chrB_fasta_file = (ref_genome_fasta, patch_db_fasta) if chrB in patch_db_entries else (patch_db_fasta, ref_genome_fasta)

        chrA_lend = coordA - pad_region_length
        chrA_rend = coordA
        if brkpt_type != 'Split':
            chrA_rend += pad_region_length
            
        if chrA_lend < 1:
            chrA_lend = 1
                    
        chrA_seq_region = extract_seq_region(chrA_fasta_file, chrA, chrA_lend, chrA_rend, orientA)

        chrB_lend = coordB
        if brkpt_type != 'Split':
            chrB_lend -= pad_region_length
        
        chrB_rend = coordB + pad_region_length
        if chrB_lend < 1:
            chrB_lend = 1

        chrB_seq_region = extract_seq_region(chrB_fasta_file, chrB, chrB_lend, chrB_rend, orientB)

        # build target sequence and gtf 
        concat_seq = chrA_seq_region
        
        print("\t".join(str(x) for x in [event_acc, "VIF-draft", "region", 1, len(concat_seq), ".", orientA, ".", "{} {}-{}".format(chrA, chrA_lend, chrA_rend)]), file=ofh_gtf)

        concat_seq += 'N' * 100 # add spacer

        print("\t".join(str(x) for x in [event_acc, "VIF-draft", "region", len(concat_seq) + 1, len(concat_seq) + len(chrB_seq_region), ".", orientB, ".", "{} {}-{}".format(chrB, chrB_lend, chrB_rend)]), file=ofh_gtf)

        concat_seq += chrB_seq_region


        print(">{}\n{}".format(event_acc, concat_seq), file=ofh_fasta)


    ofh_fasta.close()
    ofh_gtf.close()
        
    
    return

        
        

def extract_seq_region(fasta_filename, chr, lend, rend, orient):

    cmd = "samtools faidx {} {}:{}-{}".format(fasta_filename, chr, lend, rend, orient)

    if orient == '-':
        cmd += " --reverse-complement"

    result = "".join(subprocess.check_output(cmd, shell=True, encoding='utf-8').split("\n")[1:])
        
    return result


def parse_chim_events(chim_events_filename):

    event_info_dict = dict()

    csv.field_size_limit(int(1e6))
    
    with open(chim_events_filename, 'rt') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            event_num = row['entry']
            event_info_dict[event_num] = row

    return event_info_dict
    
            
    
    


if __name__=='__main__':
    main()

    
