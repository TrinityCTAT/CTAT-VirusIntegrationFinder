#!/usr/bin/env python3
# encoding: utf-8

import argparse
import glob
import os
import sys

if sys.version_info[0] < 3:
    raise Exception("This script requires Python 3")

sys.path.insert(
    0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../PyLib"])
)
from Pipeliner import Pipeliner, Command

import logging

FORMAT = "%(asctime)-15s %(levelname)s %(module)s.%(name)s.%(funcName)s at %(lineno)d :\n\t%(message)s\n"
global logger
logger = logging.getLogger()
logging.basicConfig(format=FORMAT, level=logging.DEBUG)

BASEDIR = os.path.dirname(os.path.abspath(__file__))
UTILDIR = BASEDIR
LRsplice_UTILDIR = os.path.join(BASEDIR, "LRsplice_util")


def main():
    arg_parser = argparse.ArgumentParser(
        description="Aligns to virus-integrated genome region",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    arg_parser._action_groups.pop()
    required = arg_parser.add_argument_group("required arguments")
    optional = arg_parser.add_argument_group("optional arguments")

    required.add_argument(
        "--left_fq",
        type=str,
        required=False,
        default=None,
        help="left (or single) fastq file",
    )

    required.add_argument(
        "--right_fq",
        type=str,
        required=False,
        default="",  # intentionally not None
        help="right fastq file (optional)",
    )

    required.add_argument(
        "--virus_insertions_tsv",
        type=str,
        required=True,
        help="virus insertion candidates tsv file. Only the first entry is targeted.",
    )

    required.add_argument(
        "--flank",
        type=int,
        required=True,
        help="megabases to flank the virus insertion",
    )
    
    optional.add_argument(
        "--genome_lib_dir",
        type=str,
        default=os.environ.get("CTAT_GENOME_LIB"),
        help="genome lib directory - see http://FusionFilter.github.io for details.  Uses env var CTAT_GENOME_LIB as default",
    )

    required.add_argument(
        "--viral_db_fasta", type=str, required=True, help="viral db fasta file"
    )

    optional.add_argument(
        "--viral_db_gtf",
        type=str,
        required=False,
        default=None,
        help="viral db gtf file",
    )

    optional.add_argument(
        "-O",
        "--output_dir",
        dest="output_dir",
        type=str,
        required=False,
        default="VIF.LRsplice.outdir",
        help="output directory",
    )

    optional.add_argument(
        "--out_prefix", type=str, default="vif.LRsplice", help="output filename prefix"
    )
    optional.add_argument(
        "--CPU",
        required=False,
        type=int,
        default=4,
        help="number of threads for multithreaded processes",
    )

    optional.add_argument(
        "--remove_duplicates",
        action="store_true",
        default=False,
        help="remove duplicate alignments",
    )

    args_parsed = arg_parser.parse_args()

    left_fq = os.path.abspath(args_parsed.left_fq)
    right_fq = os.path.abspath(args_parsed.right_fq) if args_parsed.right_fq else ""
    if left_fq == right_fq:
        raise ValueError("Left and right fastqs are the same.")

    output_dir = os.path.abspath(args_parsed.output_dir)
    genome_lib_dir = os.path.abspath(args_parsed.genome_lib_dir)
    viral_db_fasta = os.path.abspath(args_parsed.viral_db_fasta)
    viral_db_gtf = (
        os.path.abspath(args_parsed.viral_db_gtf) if args_parsed.viral_db_gtf else ""
    )

    virus_insertions_tsv = os.path.abspath(args_parsed.virus_insertions_tsv)
    flank_mb = args_parsed.flank

    remove_duplicates_flag = args_parsed.remove_duplicates

    output_prefix = args_parsed.out_prefix

    if remove_duplicates_flag:
        output_prefix = output_prefix + ".DupsRm"

    if not genome_lib_dir:
        sys.stderr.write(
            "Error, must set --genome_lib_dir or have env var CTAT_GENOME_LIB set"
        )
        sys.exit(1)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    os.chdir(output_dir)

    checkpoints_dir = os.path.join(
        output_dir, "__chckpts_VIF_LRsplice-rmdups-{}".format(remove_duplicates_flag),
    )
    pipeliner = Pipeliner(checkpoints_dir)

    #
    #
    # -- begin pipeline --
    #
    #


    cmd = " ".join(
        [
            os.path.join(LRsplice_UTILDIR, "prep_viral_genome_insertion_w_flank.py"),
            "--virus_insertion {}".format(virus_insertions_tsv),
            "--genome_lib_dir {}".format(genome_lib_dir),
            "--viral_db_fasta {}".format(viral_db_fasta),
            "--flank {}".format(flank_mb),
            "--output_prefix {}".format(output_prefix),
            "--output_dir {}".format(output_dir),
        ]
    )
    if viral_db_gtf:
        cmd += " --viral_db_gtf {}".format(viral_db_gtf)


    pipeliner.add_commands([Command(cmd, "prep_LRsplice")])
    
    LRsplice_fasta = os.path.join(output_dir, output_prefix + ".fasta")
    LRsplice_gtf = os.path.join(output_dir, output_prefix + ".gtf")
    LRsplice_virus_only_gtf = os.path.join(output_dir, output_prefix + ".virus-only.gtf")


    
    cmd = " ".join(
        [
            os.path.join(UTILDIR, "STAR_chimeric_patch_runner.py"),
            "--left_fq {}".format(left_fq),
            "--patch_db_fasta {}".format(LRsplice_fasta),
            "--patch_db_gtf {}".format(LRsplice_gtf),
            "--disable_chimeras",
            "--genome_lib_dir {}".format(genome_lib_dir),
            "-O " + output_dir,
        ]
    )
    
    if right_fq:
        cmd += " --right_fq {}".format(right_fq)

        


    pipeliner.add_commands(
        [Command(cmd, "star.LRsplice.rmdups-{}".format(remove_duplicates_flag))]
    )

    LRsplice_bam = os.path.join(
        output_dir, "Aligned.sortedByCoord.out.bam"
    )
    
    if remove_duplicates_flag:
        ## remove duplicate alignments
        pipeliner.add_commands(
            [Command("samtools index " + LRsplice_bam, "index_star_bam",)]
        )

        LRsplice_bam_dups_removed = os.path.join(
            outdir, "Aligned.sortedByCoord.out.dups_removed.bam"
        )
        cmd = " ".join(
            [
                os.path.join(UTILDIR, "bam_mark_duplicates.py"),
                "-i {}".format(LRsplice_bam),
                "-o {}".format(LRsplice_bam_dups_removed),
                "-r",
            ]
        )
        pipeliner.add_commands([Command(cmd, "LRsplice_bam_rmdups")])
        pipeliner.add_commands(
            [
                Command(
                    "samtools index " + LRsplice_bam_dups_removed,
                    "index_LRsplice_rmdups_bam",
                )
            ]
        )

        LRsplice_bam = LRsplice_bam_dups_removed


    ## Run pipeline
    pipeliner.run()

    sys.exit(0)


if __name__ == "__main__":
    main()
