#!/usr/bin/env python3
# encoding: utf-8

import os, re, sys
import argparse
import subprocess
import pysam
import csv

if sys.version_info[0] < 3:
    raise Exception("This script requires Python 3")

sys.path.insert(
    0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../PyLib"])
)

import logging

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)


def main():

    arg_parser = argparse.ArgumentParser(
        description="models genome/viral insertion for examining long-range virus/host splicing",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    arg_parser.add_argument(
        "--genome_lib_dir",
        type=str,
        default=os.environ.get("CTAT_GENOME_LIB"),
        help="genome lib directory - see http://FusionFilter.github.io for details.  Uses env var CTAT_GENOME_LIB as default",
    )

    arg_parser.add_argument(
        "--viral_db_fasta", type=str, required=True, help="viral genome database"
    )

    arg_parser.add_argument(
        "--viral_db_gtf",
        type=str,
        required=False,
        default=None,
        help="viral genome annotations in gtf format",
    )

    arg_parser.add_argument(
        "--output_dir", type=str, required=True, help="output directory"
    )

    arg_parser.add_argument(
        "--output_prefix",
        type=str,
        required=True,
        help="output prefix for fasta and gtf files",
    )

    arg_parser.add_argument(
        "--virus_insertion",
        type=str,
        required=True,
        help="virus insertion candidates tsv file (only the first entry is processed)",
    )

    arg_parser.add_argument(
        "--flank",
        type=int,
        required=True,
        help="length (in MB) of genome regions flanking the virus insertion",
    )

    args_parsed = arg_parser.parse_args()

    genome_lib_dir = args_parsed.genome_lib_dir
    viral_db_fasta = args_parsed.viral_db_fasta
    viral_db_gtf = args_parsed.viral_db_gtf
    output_dir = args_parsed.output_dir
    output_prefix = args_parsed.output_prefix
    flank_mb = args_parsed.flank
    virus_insertion_tsv = args_parsed.virus_insertion

    if not genome_lib_dir:
        logger.error("Error, --genome_lib_dir must be specified")
        sys.exit(1)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    insertion_event = parse_insertion_event(virus_insertion_tsv)

    genome_left_struct, virus_struct, genome_right_struct = model_virus_insertion(
        insertion_event, flank_mb
    )

    build_virus_insertion_sequence(
        genome_left_struct,
        virus_struct,
        genome_right_struct,
        genome_lib_dir,
        viral_db_fasta,
        viral_db_gtf,
        output_dir,
        output_prefix,
    )

    sys.exit(0)


def build_virus_insertion_sequence(
    genome_left_struct,
    virus_struct,
    genome_right_struct,
    genome_lib_dir,
    viral_db_fasta,
    viral_db_gtf,
    output_dir,
    output_prefix,
):

    ref_genome_fasta = os.path.join(genome_lib_dir, "ref_genome.fa")
    ref_annot_gtf = os.path.join(genome_lib_dir, "ref_annot.gtf")

    out_fasta_filename = os.path.join(output_dir, output_prefix + ".fasta")
    out_gtf_filename = os.path.join(output_dir, output_prefix + ".gtf")
    out_gtf_virus_only_filename = os.path.join(output_dir, output_prefix + ".virus-only.gtf")

    ofh_fasta = open(out_fasta_filename, "wt")
    ofh_gtf = open(out_gtf_filename, "wt")
    ofh_virus_gtf = open(out_gtf_virus_only_filename, "wt")
    
    left_host_gtf, right_host_gtf = extract_host_genome_annotations(
        genome_left_struct, genome_right_struct, ref_annot_gtf
    )
    
    left_genomic_region_seq = extract_seq_region(
        ref_genome_fasta,
        genome_left_struct["chrom"],
        genome_left_struct["lend"],
        genome_left_struct["rend"],
        "+",
    )
    right_genomic_region_seq = extract_seq_region(
        ref_genome_fasta,
        genome_right_struct["chrom"],
        genome_right_struct["lend"],
        genome_right_struct["rend"],
        "+",
    )

    virus_genome_seq = extract_viral_genome(
        viral_db_fasta, virus_struct["virus"], virus_struct["orient"]
    )

    virus_genome_gtf = (
        extract_virus_genome_annotations(viral_db_gtf, virus_struct["virus"])
        if viral_db_gtf
        else ""
    )

    ## build virus insertion model, output seq and gtf

    # left host
    left_genome_start = genome_left_struct["lend"]
    concat_seq = left_genomic_region_seq
    for gtf_line in left_host_gtf:
        vals = gtf_line.split("\t")
        vals[0] = "imodel"
        vals[3] = str(int(vals[3]) - left_genome_start + 1)
        vals[4] = str(int(vals[4]) - left_genome_start + 1)

        print("\t".join(vals), file=ofh_gtf)

    # virus
    pseudo_start = len(concat_seq)
    virus_genome_seq_len = len(virus_genome_seq)

    ## add entry for just the length of the virus for easy identification:
    virus_region = ["imodel", "virus", "region",
                     str(pseudo_start + 1), str(pseudo_start + virus_genome_seq_len),
                     ".", virus_struct["orient"], ".",
                     "gene_id \"{}\"".format(virus_struct["virus"]) ]
    print("\t".join(virus_region), file=ofh_gtf)
    print("\t".join(virus_region), file=ofh_virus_gtf)
    

    for gtf_line in virus_genome_gtf:
        vals = gtf_line.split("\t")
        vals[0] = "imodel"
        if virus_struct["orient"] == "-":
            # revcomp the virus coords
            virus_lend, virus_rend = int(vals[3]), int(vals[4])
            virus_lend, virus_rend = sorted(
                [
                    virus_genome_seq_len - virus_lend + 1,
                    virus_genome_seq_len - virus_rend + 1,
                ]
            )
            vals[3] = str(virus_lend)
            vals[4] = str(virus_rend)
            vals[6] = "-" if vals[6] == "+" else "+"

        vals[3] = str(pseudo_start + int(vals[3]))
        vals[4] = str(pseudo_start + int(vals[4]))

        print("\t".join(vals), file=ofh_gtf)
        print("\t".join(vals), file=ofh_virus_gtf)
        
    # add virus seq
    concat_seq += virus_genome_seq

    # right host
    right_genome_start = genome_right_struct["lend"]
    pseudo_start = len(concat_seq)
    for gtf_line in right_host_gtf:
        vals = gtf_line.split("\t")
        vals[0] = "imodel"
        vals[3] = str(int(vals[3]) - right_genome_start + pseudo_start +1)
        vals[4] = str(int(vals[4]) - right_genome_start + pseudo_start +1)

        print("\t".join(vals), file=ofh_gtf)

    concat_seq += right_genomic_region_seq

    print(">{}\n{}".format("imodel", concat_seq), file=ofh_fasta)

    ofh_fasta.close()
    ofh_gtf.close()
    ofh_virus_gtf.close()
    
    return


def extract_seq_region(fasta_filename, chr, lend, rend, orient):

    cmd = "samtools faidx {} {}:{}-{}".format(fasta_filename, chr, lend, rend, orient)

    if orient == "-":
        cmd += " --reverse-complement"

    result = "".join(
        subprocess.check_output(cmd, shell=True, encoding="utf-8").split("\n")[1:]
    )

    return result


def extract_viral_genome(fasta_filename, virus, orient):

    cmd = "samtools faidx {} {}".format(fasta_filename, virus, orient)

    if orient == "-":
        cmd += " --reverse-complement"

    result = "".join(
        subprocess.check_output(cmd, shell=True, encoding="utf-8").split("\n")[1:]
    )

    return result


def extract_host_genome_annotations(
    genome_left_struct, genome_right_struct, ref_annot_gtf
):

    left_genome_annots = list()
    right_genome_annots = list()

    with open(ref_annot_gtf) as fh:
        for line in fh:
            line = line.rstrip()
            if line[0] == "#":
                continue
            vals = line.split("\t")
            if len(vals) < 8:
                continue
            
            newline = feature_in_range(vals, genome_left_struct)
            if newline:
                left_genome_annots.append(newline)

            newline = feature_in_range(vals, genome_right_struct)
            if newline:
                right_genome_annots.append(newline)

    return left_genome_annots, right_genome_annots


def extract_virus_genome_annotations(viral_db_gtf, virus_acc):

    virus_annots = list()

    with open(viral_db_gtf) as fh:
        for line in fh:
            
            line = line.rstrip()
            vals = line.split("\t")

            if vals[0] == virus_acc:
                virus_annots.append(line)

    if not virus_annots:
        logger.warn(
            "No virus annotations found for {} in {}".format(virus_acc, viral_db_gtf)
        )

    return virus_annots


def feature_in_range(vals, genome_struct):

    vals = list(vals)  # make copy

    lend = int(vals[3])
    rend = int(vals[4])

    newline = None

    if (vals[0] == genome_struct['chrom'] and
        lend < genome_struct["rend"] and
        rend > genome_struct["lend"]):
        
        # got overlap

        if lend < genome_struct["lend"]:
            vals[3] = str(genome_struct["lend"])

        if rend > genome_struct["rend"]:
            vals[4] = str(genome_struct["rend"])

        newline = "\t".join(vals)

    return newline


def parse_insertion_event(insertion_event_filename):

    event_info_dict = dict()

    csv.field_size_limit(int(1e6))

    with open(insertion_event_filename, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            return row
            # event_num = row["entry"]
            # event_info_dict[event_num] = row

    raise RuntimeError(
        "Error, no insertion candidates identified in {}".format(chim_events_filename)
    )


def model_virus_insertion(insertion_event, flank_mb):

    chrA = insertion_event["chrA"]
    coordA = int(insertion_event["coordA"])
    orientA = insertion_event["orientA"]

    chrB = insertion_event["chrB"]
    coordB = int(insertion_event["coordB"])
    orientB = insertion_event["orientB"]

    flank_bp = int(int(flank_mb) * 1e6)

    if re.match("chr", chrA):
        assert orientA == "+", "Error, host genome chr is not in + orient as expected"
        genome_left_struct = {
            "chrom": chrA,
            "lend": max(1, coordA - flank_bp),
            "rend": coordA,
        }

        virus_struct = {"virus": chrB, "orient": orientB}

        genome_right_struct = {
            "chrom": chrA,
            "lend": coordA + 1,
            "rend": coordA + 1 + flank_bp,
        }

    else:
        assert re.match(
            "chr", chrB
        ), "Error, neither chrom of insertion candidate begins with 'chr' "

        assert orientB == "+", "Error, host genome chr is not in + orient as expected"

        genome_left_struct = {
            "chrom": chrB,
            "lend": max(1, coordB - 1 - flank_bp),
            "rend": coordB - 1,
        }

        virus_struct = {"virus": chrA, "orient": orientA}

        genome_right_struct = {"chrom": chrB, "lend": coordB, "rend": coordB + flank_bp}

    return genome_left_struct, virus_struct, genome_right_struct


if __name__ == "__main__":
    main()
