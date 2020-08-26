#!/usr/bin/env python3
# encoding: utf-8

import os, re, sys
import argparse
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

UTILDIR = os.path.dirname(__file__)


def main():
    
    arg_parser = argparse.ArgumentParser(description="writes star commands for insertion event evaluations",
                                         formatter_class=argparse.RawTextHelpFormatter
                                         )
    arg_parser.add_argument("--workdir_base", type=str, required=True,
                            help="base directory for work")
    
    arg_parser.add_argument("--output_file", "-O", dest="output_filename", type=str, required=True,
                            help="output file to write commands to")
    
        
    args_parsed = arg_parser.parse_args()


    workdir_base = os.path.abspath(args_parsed.workdir_base)
    output_filename = os.path.abspath(args_parsed.output_filename)

    ofh = open(output_filename, 'wt')    

    chim_targets_file = os.path.join(workdir_base, "chim_events_for_eval.tsv")

    with open(chim_targets_file, 'rt') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            workdir = row['workdir']
            
            target_genome_fa = os.path.join(workdir, "target.fasta")
            left_fq_filename = os.path.join(workdir, "reads_R1.fastq")
            right_fq_filename = os.path.join(workdir, "reads_R2.fastq")

            cmd = " ".join([ os.path.join(UTILDIR, "STAR_nonchimeric_patchless_runner.py"),
                             "--genome_fa {}".format(target_genome_fa),
                             "--left_fq {}".format(left_fq_filename),
                             "--CPU 1"])

            if os.path.exists(right_fq_filename):
                cmd += "--right_fq {}".format(right_fq_filename)

            print(cmd, file=ofh)

    ofh.close()

    logger.info("-wrote star eval cmds {}".format(output_filename))
    
    sys.exit(0)


if __name__=='__main__':
    main()
