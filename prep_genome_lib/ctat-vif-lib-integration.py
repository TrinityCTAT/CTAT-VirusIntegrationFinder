#!/usr/bin/env python3

import sys, os, re
import subprocess
import argparse
import logging
import pandas as pd

sys.path.insert(
    0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "PyLib"])
    )
from Pipeliner import Pipeliner, Command

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

    checkpoints_dir = f"{VIF_dir}/__checkpts.dir"
    pipeliner = Pipeliner(checkpoints_dir)
    
    # copy the virus db to the VIF dir and index it.
    installed_virus_db = os.path.join(VIF_dir, "virus_db.fasta")
    logger.info("-installing virus db")
    pipeliner.add_commands([Command(f"cp {virus_db_fasta} {installed_virus_db}", "cp_virus_to_VIF.ok")])
    pipeliner.add_commands([Command(f"samtools faidx {installed_virus_db}", "faidx_virusdb.ok")])
    pipeliner.run()

    # mask homologous regions from ref genome:
    logger.info("-masking virus homologous regions from ref genome")
    cmdstr = f"makeblastdb -in {ref_genome_fasta} -dbtype nucl"
    pipeliner.add_commands([Command(cmdstr, "makeblastableref.ok")])
    pipeliner.run()

    blastn_outfile = f"{VIF_dir}/genome_virus_blastn.outfmt6"
    cmdstr = f"blastn -query {virus_db_fasta} -db {ref_genome_fasta} -outfmt 6 -evalue 1e-10 -num_threads {num_threads} > {blastn_outfile}"
    pipeliner.add_commands([Command(cmdstr, "blastnvirustogenome.ok")])
    pipeliner.run()

    ## generate masked genome.
    # make regions file
    if sum(1 for l in open(blastn_outfile))>0:
        df = pd.read_csv(blastn_outfile, sep="\t", header=None, usecols=[1,2,8,9])
        df.columns = ["chrom", "perid", "end5", "end3"]
        df = df[ df['perid'] >= 90 ]  # filter, require min 90% identity
        df['lend'] = df.apply(lambda x: min(x['end5'], x['end3']), axis=1)
        df['rend'] = df.apply(lambda x: max(x['end5'], x['end3']), axis=1)
        df['lend'] = df['lend']-1 
        df = df[ ['chrom', 'lend', 'rend'] ]

    else: 
        df = pd.DataFrame(columns=['chrom', 'lend', 'rend'])
        
    bed_filename = blastn_outfile + ".bed"
    df.to_csv(bed_filename, sep="\t", header=False, index=False)
    
    # do the masking using bedtools
    virus_masked_genome = f"{VIF_dir}/ref_genome.virus_masked.fa"
    cmdstr = f"maskFastaFromBed -fi {ref_genome_fasta} -fo {virus_masked_genome} -bed {bed_filename}"
    pipeliner.add_commands([Command(cmdstr, "virusmaskingrefgenome.ok")])
    pipeliner.run()

    # before masking: 151122963
    # after masking:  152498464
    # so, 1,375,501 new bases masked.
    
    
    # create new fasta file including viruses and human genome together:
    combined_genomes_fa = os.path.join(VIF_dir, "hg_plus_viraldb.fasta")
    logger.info(f"-combining {virus_masked_genome} and {virus_db_fasta} into {combined_genomes_fa}")
    cmdstr = f"cat {virus_masked_genome} {virus_db_fasta} > {combined_genomes_fa}"
    pipeliner.add_commands([Command(cmdstr, "combineMaskedGenomeWithViruses.ok")])
        
    cmdstr = f"samtools faidx {combined_genomes_fa}"
    pipeliner.add_commands([Command(cmdstr, "combinedgenomes.faidx.ok")])
    pipeliner.run()

        
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


    pipeliner.add_commands([Command(star_index_cmd, "buildCombinedIndex.ok")])
    pipeliner.run()
    
    logger.info("-Done.")

    sys.exit(0)



if __name__=='__main__':
    main()

