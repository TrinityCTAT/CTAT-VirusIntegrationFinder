#!/usr/bin/env python3
# encoding: utf-8

import os, re, sys
import argparse
import subprocess

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
    arg_parser._action_groups.pop()
    required = arg_parser.add_argument_group('required arguments')
    optional = arg_parser.add_argument_group('optional arguments')
    
    required.add_argument("--left_fq", type=str, required=False, default=None,
                          help="left (or single) fastq file")
    
    required.add_argument("--right_fq", type=str, required=False, default="", #intentionally not None
                          help="right fastq file (optional)")
    
    optional.add_argument("--genome_lib_dir", type=str, default=os.environ.get('CTAT_GENOME_LIB'),
                          help="genome lib directory - see http://FusionFilter.github.io for details.  Uses env var CTAT_GENOME_LIB as default")
    
    optional.add_argument("-O", "--output_dir", type=str, required=False, default="VIF.outdir", help="output directory")
    optional.add_argument("--CPU", dest="CPU", required=False, type=int, default=4,
                          help="number of threads for multithreaded processes")

    required.add_argument("--patch_db_fasta", type=str, required=True,
                             help="database of additional genome targets to add")

    optional.add_argument("--patch_db_gtf", type=str, required=False,
                          help="gtf annotations for patch_db_fasta entries")

    optional.add_argument("--max_mate_dist", type=int, required=False, default=100000,
                          help="max distance between mates and max intron length allowed")
    
    args_parsed = arg_parser.parse_args()
    
    
    left_fq = os.path.abspath(args_parsed.left_fq)
    right_fq = os.path.abspath(args_parsed.right_fq) if args_parsed.right_fq else ""
    output_dir = os.path.abspath(args_parsed.output_dir)
    genome_lib_dir = args_parsed.genome_lib_dir
    patch_db_fasta = os.path.abspath(args_parsed.patch_db_fasta)
    patch_db_gtf = os.path.abspath(args_parsed.patch_db_gtf) if args_parsed.patch_db_gtf else ""
    
    if not genome_lib_dir:
        sys.stderr.write("Error, must set --genome_lib_dir or have env var CTAT_GENOME_LIB set")
        sys.exit(1)
        
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    os.chdir(output_dir)

    checkpoints_dir = os.path.join(output_dir, "__chckpts_star_chimeric")
    pipeliner = Pipeliner(checkpoints_dir)


    readFilesIn = left_fq
    if right_fq:
        readFilesIn += " " + right_fq
    
    cmd = " ".join(["STAR",
                    "--runThreadN {}".format(args_parsed.CPU),
                    "--genomeDir {}".format(os.path.join(genome_lib_dir, "ref_genome.fa.star.idx")),
                    "--outSAMtype BAM SortedByCoordinate",
                    "--twopassMode Basic",
                    "--alignSJDBoverhangMin 10",
                    "--genomeSuffixLengthMax 10000",
                    "--limitBAMsortRAM 47271261705", 
                    "--alignInsertionFlush Right",
                    "--alignMatesGapMax {}".format(args_parsed.max_mate_dist),
                    "--alignIntronMax {}".format(args_parsed.max_mate_dist),
                    "--readFilesIn {}".format(readFilesIn),
                    "--chimJunctionOverhangMin 12",
                    "--chimSegmentMin 12",
                    "--chimSegmentReadGapMax 3",
                    "--genomeFastaFiles {}".format(patch_db_fasta),
                    "--outSAMfilter KeepAllAddedReferences",
                    "--alignSJstitchMismatchNmax 5 -1 5 5",
                    "--scoreGapNoncan -6" ])


    if patch_db_gtf:
        cmd += " --sjdbGTFfile {} --sjdbOverhang 150 "

    if left_fq[-3:] == ".gz":
        cmd += " --readFilesCommand 'gunzip -c' "

    
    pipeliner.add_commands([Command(cmd, "star_chimeric")])
    
    
    pipeliner.run()

    sys.exit(0)


if __name__=='__main__':
    main()
