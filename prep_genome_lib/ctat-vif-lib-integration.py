#!/usr/bin/env python3

import sys, os, re
import subprocess
import argparse
import logging


import logging
FORMAT = "%(asctime)-15s: %(levelname)s %(module)s.%(name)s.%(funcName)s %(message)s"
logger = logging.getLogger(__file__)
logging.basicConfig(stream=sys.stderr, format=FORMAT, level=logging.INFO)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--virus_db",
                        type=str, required=True,
                        help="virus database (multi-fasta file)")

    parser.add_argument("--genome_lib_dir", dest="genome_lib_dir",
                        type=str, required=True,
                        help="CTAT genome lib build dir")

    parser.add_argument("--CPU",
                        type=int,
                        default=4,
                        help="num threads to use for building STAR index")

    args=parser.parse_args()
    
    genome_lib_dir = args.genome_lib_dir
    virus_db_fasta = args.virus_db
    num_threads = args.CPU
    
    ref_genome_fasta = os.path.join(genome_lib_dir, "ref_genome.fa")
    ref_gtf = os.path.join(genome_lib_dir, "ref_annot.gtf")

    if not os.path.exists(ref_genome_fasta):
        exit("Error, not finding genome at: {}".format(ref_genome_fasta))

    if not os.path.exists(ref_gtf):
        exit("Error, not finding ref annotation at: {}".format(ref_gtf))


    VIF_dir = os.path.join(genome_lib_dir, "VIF")
    if not os.path.exists(VIF_dir):
        os.makedirs(VIF_dir)


    def run_cmd(cmd):
        logger.info(f"CMD: {cmd}")
        subprocess.check_call(cmd, shell=True)        


    # copy the virus db to the VIF dir and index it.
    installed_virus_db = os.path.join(VIF_dir, "virus_db.fasta")
    logger.info("-installing virus db")
    run_cmd(f"cp {virus_db} {installed_virus_db}")
    run_cmd(f"samtools faidx {installed_virus_db}")
    

    # create new fasta file including viruses and human genome together:
    logger.info(f"-combining {ref_genome_fasta} and {virus_db} into {combined_genomes_fa}")
    combined_genomes_fa = os.path.join(VIF_dir, "hg_plus_viraldb.fasta")
    run_cmd(f"cat {ref_genome_fasta} {virus_db} > {combined_genomes_fa}")
    run_cmd(f"samtools faidx {combined_genomes_fa}")
    
    # build star index
    logger.info(f"-building star index for {combined_genomes_fa}")
    star_index_cmd = " ".join(["STAR",
                               "--runMode genomeGenerate",
                               "--runThreadN {}".format(num_threads),
                               "--genomeDir {}.star.idx".format(combined_genomes_fa),
                               "--genomeFastaFiles {}".format(combined_genomes_fa),
                               "--limitGenomeGenerateRAM 40419136213",
                               "--sjdbGTFfile {}".format(ref_gtf),
                               "--sjdbOverhang 100"])
    run_cmd(star_index_cmd)
    

    logger.info("-Done.")

    sys.exit(0)



if __name__=='__main__':
    main()

