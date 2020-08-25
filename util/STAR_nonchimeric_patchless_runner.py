#!/usr/bin/env python3
# encoding: utf-8

import os, re, sys
import argparse
import subprocess
import math

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
    
    required.add_argument("--left_fq", type=str, required=True, default=None,
                          help="left (or single) fastq file")
    
    required.add_argument("--right_fq", type=str, required=False, default="", #intentionally not None
                          help="right fastq file (optional)")
    
    optional.add_argument("--genome_fa", type=str, required=True,
                          help="target genome fasta file")
    
    optional.add_argument("-O", "--output_dir", type=str, required=False, default="star.outdir", help="output directory")
    optional.add_argument("--CPU", dest="CPU", required=False, type=int, default=4,
                          help="number of threads for multithreaded processes")

    
    optional.add_argument("--max_mate_dist", type=int, required=False, default=100000,
                          help="max distance between mates and max intron length allowed")
    
    args_parsed = arg_parser.parse_args()
    
    
    left_fq = os.path.abspath(args_parsed.left_fq)
    right_fq = os.path.abspath(args_parsed.right_fq) if args_parsed.right_fq else ""
    output_dir = os.path.abspath(args_parsed.output_dir)
    genome_fa = os.path.abspath(args_parsed.genome_fa)
            
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    os.chdir(output_dir)

    checkpoints_dir = os.path.join(output_dir, "__chckpts_star_nonchimeric")
    pipeliner = Pipeliner(checkpoints_dir)



    genome_size = os.path.getsize(genome_fa)
    min_ram = 1024**3
    estimated_ram = genome_size * 15
    estimated_ram = max(min_ram, estimated_ram)

    genomeSAindexNbases = int(math.log2(genome_size) / 2)


    starGenomeIdxDir = "target.fasta.star.idx"
    if not os.path.exists(starGenomeIdxDir):
        os.makedirs(starGenomeIdxDir)
    
    ## genome generate:
    cmd = " ".join(["STAR",
                    "--runThreadN {}".format(args_parsed.CPU),
                    "--runMode genomeGenerate",
                    "--genomeDir {}".format(starGenomeIdxDir),
                    "--genomeFastaFiles {}".format(genome_fa),
                    "--genomeSAindexNbases {}".format(genomeSAindexNbases),
                    "--limitGenomeGenerateRAM {}".format(estimated_ram)])

    pipeliner.add_commands([Command(cmd, "star_nonchim_pless_genomeGenerate")])
    
    

    readFilesIn = left_fq
    if right_fq:
        readFilesIn += " " + right_fq
    
    cmd = " ".join(["STAR",
                    "--runThreadN {}".format(args_parsed.CPU),
                    "--genomeDir {}".format(starGenomeIdxDir),
                    "--outSAMtype BAM SortedByCoordinate",
                    "--twopassMode Basic",
                    "--alignSJDBoverhangMin 10",
                    "--genomeSuffixLengthMax 10000",
                    "--limitBAMsortRAM 47271261705", 
                    "--alignInsertionFlush Right",
                    "--alignMatesGapMax {}".format(args_parsed.max_mate_dist),
                    "--alignIntronMax {}".format(args_parsed.max_mate_dist),
                    "--readFilesIn {}".format(readFilesIn),
                    "--alignSJstitchMismatchNmax 5 -1 5 5",
                    "--scoreGapNoncan -6" ])


    if left_fq[-3:] == ".gz":
        cmd += " --readFilesCommand 'gunzip -c' "

    
    pipeliner.add_commands([Command(cmd, "star_nonchim_pless_align")])
    
    
    pipeliner.run()

    sys.exit(0)


if __name__=='__main__':
    main()
