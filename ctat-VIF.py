#!/usr/bin/env python3
# encoding: utf-8

import os, re, sys
import argparse
import subprocess

if sys.version_info[0] < 3:
    raise Exception("This script requires Python 3")

sys.path.insert(0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "PyLib"]))
from Pipeliner import Pipeliner, Command
    

import logging
FORMAT = "%(asctime)-15s %(levelname)s %(module)s.%(name)s.%(funcName)s at %(lineno)d :\n\t%(message)s\n"
global logger
logger = logging.getLogger()
logging.basicConfig(filename='ctat-VIF.log', format=FORMAT, filemode='w', level=logging.DEBUG)
# add a new Handler to print all INFO and above messages to stdout
ch = logging.StreamHandler(sys.stdout)
ch.setLevel(logging.INFO)
logger.addHandler(ch)


VERSION = "__BLEEDING_EDGE__"


BASEDIR = os.path.dirname(os.path.abspath(__file__))
UTILDIR = os.sep.join([BASEDIR, "util"])




def main():

    arg_parser = argparse.ArgumentParser(
        
    description = "Finds virus integration sites in the genome",
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
    
    required.add_argument("--viral_db_fasta", type=str, required=True,
                          help="viral db fasta file")
    optional.add_argument("--viral_db_gtf", type=str, required=False, default=None,
                           help="viral db gtf file")

        
    optional.add_argument("-O", "--output_dir", dest="output_dir", type=str, required=False, default="VIF.outdir", help="output directory")
    
    optional.add_argument("--out_prefix", type=str, default='vif', help="output filename prefix")
    optional.add_argument("--CPU", required=False, type=int, default=4,
                          help="number of threads for multithreaded processes")

    optional.add_argument("--remove_duplicates", action='store_true', default=False, help="remove duplicate alignments")
    
    args_parsed = arg_parser.parse_args()
    

    left_fq = os.path.abspath(args_parsed.left_fq)
    right_fq = os.path.abspath(args_parsed.right_fq) if args_parsed.right_fq else ""
    output_dir = os.path.abspath(args_parsed.output_dir)
    genome_lib_dir = os.path.abspath(args_parsed.genome_lib_dir)
    viral_db_fasta = os.path.abspath(args_parsed.viral_db_fasta)
    viral_db_gtf = os.path.abspath(args_parsed.viral_db_gtf) if args_parsed.viral_db_gtf else ""

    remove_duplicates_flag = args_parsed.remove_duplicates

    output_prefix = args_parsed.out_prefix

    if remove_duplicates_flag:
        output_prefix = output_prefix + ".DupsRm"
    
    if not genome_lib_dir:
        sys.stderr.write("Error, must set --genome_lib_dir or have env var CTAT_GENOME_LIB set")
        sys.exit(1)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    os.chdir(output_dir)

    checkpoints_dir = os.path.join(output_dir, "__chckpts_VIF")
    pipeliner = Pipeliner(checkpoints_dir)
    
    ## run STAR using all viruses, select viruses
    cmd = " ".join([ os.path.join(UTILDIR, "STAR_chimeric_patch_runner.py"),
                     "--left_fq {}".format(left_fq),
                     "--patch_db_fasta {}".format(viral_db_fasta),
                     "--genome_lib_dir {}".format(genome_lib_dir),
                     "-O VIF_starChim_init"])
    
    if right_fq:
        cmd += " --right_fq {}".format(right_fq)
    
    pipeliner.add_commands([Command(cmd, "star_chimeric_initial")])

    local_output_prefix = "prelim.vif.dups_rm={}".format(remove_duplicates_flag)
    
    ## run virus integration site analysis, genereate report.
    cmd = " ".join([ os.path.join(UTILDIR, "chimJ_to_virus_insertion_candidate_sites.py"),
                     "--chimJ {}".format("VIF_starChim_init/Chimeric.out.junction"),
                     "--patch_db_fasta {}".format(viral_db_fasta),
                     "--output_prefix {}".format(local_output_prefix) ])

    if remove_duplicates_flag:
        cmd += " --remove_duplicates"
    
    pipeliner.add_commands([Command(cmd, "chimJ_to_insertion_candidates.dups_removed={}".format(remove_duplicates_flag))])
    
    
    ## extract targets for review
    cmd = " ".join([ os.path.join(UTILDIR, "extract_chimeric_genomic_targets.py"),
                     "--genome_lib_dir {}".format(genome_lib_dir),
                     "--patch_db_fasta {}".format(viral_db_fasta),
                     "--output_prefix {}".format("extracted.pad1k"),
                     "--chim_events {}.abridged.tsv".format(local_output_prefix),
                     "--pad_region_length 1000"])
    
    pipeliner.add_commands([Command(cmd, "extract_draft_candidates_fasta_n_gtf")])
    

    ## run STAR chimeric patch runner again, but now use candidates as patch
    
    cmd = " ".join([ os.path.join(UTILDIR, "STAR_chimeric_patch_runner.py"),
                     "--left_fq {}".format(left_fq),
                     "--patch_db_fasta {}".format("extracted.pad1k.fasta"),
                     "--disable_chimeras",
                     "--genome_lib_dir {}".format(genome_lib_dir),
                     "-O VIF_starChim_prelim_candidates"])
    
    if right_fq:
        cmd += " --right_fq {}".format(right_fq)

    pipeliner.add_commands([Command(cmd, "star_chimeric_prelim_candidates")]) 


    chimeric_bam = "VIF_starChim_prelim_candidates/Aligned.sortedByCoord.out.bam"
    
    if remove_duplicates_flag:
        ## remove duplicate alignments
        pipeliner.add_commands([Command("samtools index VIF_starChim_prelim_candidates/Aligned.sortedByCoord.out.bam", "index_star_bam")]) 

        chimeric_bam_dups_removed = "VIF_starChim_prelim_candidates/Aligned.sortedByCoord.out.dups_removed.bam"
        cmd = " ".join([os.path.join(UTILDIR, "bam_mark_duplicates.py"),
                        "-i {}".format(chimeric_bam),
                        "-o {}".format(chimeric_bam_dups_removed),
                        "-r" ] )
        pipeliner.add_commands([Command(cmd, "star_chimeric_prelim_rmdups")])
        pipeliner.add_commands([Command("samtools index VIF_starChim_prelim_candidates/Aligned.sortedByCoord.out.dups_removed.bam", "index_star_rmdups_bam")]) 

        chimeric_bam = chimeric_bam_dups_removed

    
    ## score alignments.
    scored_alignments_file = "VIF_starChim_prelim_candidates/VIF.evidence_counts.rmdups={}.tsv".format(remove_duplicates_flag)
    cmd = " ".join([os.path.join(UTILDIR, "chimeric_contig_evidence_analyzer.py"),
                    "--patch_db_bam {}".format(chimeric_bam),
                    "--patch_db_gtf extracted.pad1k.gtf",
                    "> {}".format(scored_alignments_file) ])
    pipeliner.add_commands([Command(cmd, "chim_contig_evidence_counts.rmdups={}".format(remove_duplicates_flag))])
    
    
    ## Run pipeline
    pipeliner.run()

    
    sys.exit(0)


if __name__=='__main__':
    main()
