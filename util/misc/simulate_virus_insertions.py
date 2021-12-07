#!/usr/bin/env python
# encoding: utf-8


'''
My strawman for simulating data is like so:

take our viral database and human genome as input

for (1..1000)
     select a random 1k sequence from the human genome
     select a random 1k sequence from a randomly selected virus in the database

     fuse them together, choosing the junction ends randomly and setting the orientations randomly

     simulate reads at 50x coverage using wgsim to generate fastqs
     
     Then merge all the mini fastqs together, join with some rna-seq sample such as our K562 favorite
     
     Run ctat-vif on the merged fastq
     Then, see how many of the 1k insertion sites we identify correctly.
     This involves both genome insertion site position and virus identity


# notes:
For pulling a sequence from the genome, you could use:
samtools faidx ref_genome.fa  ${chromosome}:${rand_start}-${rand_end} (edited)
and similarly for the virus, given whatever virus accession you choose
The 1..1000 is just a ballpark.   Maybe we instead just ensure that each virus in the viral dataset is represented by at least some number of insertions (ie 5 or 10  insertions) and see how many total targets that gives us.
Just remember your bookkeeping, so you know where your simulated insertion sites are and what viruses they represent.  You could also rename the simulated reads so they include this info.




# additional notes on running

Useful to chunk up the set of viruses and generate simulated insertions and fastq files for each chunk:

How I do this for running in Terra:

##  make chunks

fasta_file_chunker.pl ../Virus_db_Dec072021.nonUnq1kbMsk.filt90.fasta 100 Virus_db_Dec072021

## sim chunks

for file in *.fasta; do ~/GITHUB/CTAT_VIF/CTAT-VirusIntegrationFinder/util/misc/simulate_virus_insertions.py --virus_db $file --ins_per_virus 10 --out_prefix $file 2>&1 | tee $file.log; done

 gzip *fq

 ##  upload fastqs to terra:

 for file in *gz; do gsutil cp $file gs://fc-f5038732-f718-4441-8312-23d4e01e0882/Simulated_Data_Dec072021/$file; done

 ##  prep terra samples file

 perl -e 'print join("\n", <*fasta>);' | perl -lane '$fa = $_; $sample_id = $fa; $sample_id =~ s/\.fasta//; $gs = "gs://fc-f5038732-f718-4441-8312-23d4e01e0882/Simulated_Data_Dec072021"; print join("\t", $sample_id, "$gs/$fa.insertion_seqs.left_fq.gz", "$gs/$fa.insertion_seqs.right_fq.gz");' | tee samples.tsv

 
then manually add the samples file header as:
entity:sample_id(tab)left_fq(tab)right_fq


'''



import os, sys, re
import logging
import argparse
import subprocess
import random
import math

if sys.version_info[0] != 3:
    print("This script requires Python 3")
    exit(1)
            

logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)



def main():

    parser = argparse.ArgumentParser(description="simulate virus insertions and corresponding paired fastqs",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)


    # genome lib dir
    # virus db
    # count random targets per virus
    # read length
    # depth of coverage
    
    parser.add_argument("--genome_lib_dir", type=str, default=os.environ.get("CTAT_GENOME_LIB", None), help="ctat genome lib dir")
    parser.add_argument("--virus_db", type=str, required=True, help="virus database")
    parser.add_argument("--ins_per_virus", type=int, required=True, help="number of insertions per virus")
    parser.add_argument("--out_prefix", type=str, default="sim_vir_ins", help="output prefix for simulation fastqs and target sequences")
    parser.add_argument("--debug", action='store_true', default=False, help="debug mode")

    # can make these params sometime
    seq_extraction_len = 1000 # take 1kb from virus or ref genome for simulating reads
    depth_of_coverage = 50
    read_length = 48



    args = parser.parse_args()
    #args, unknown_args = parser.parse_known_args()


    if args.debug:
        logger.setLevel(logging.DEBUG)      
    


    genome_lib_dir = args.genome_lib_dir
    virus_db = args.virus_db
    ins_per_virus = int(args.ins_per_virus)
    out_prefix = args.out_prefix

    if genome_lib_dir is None:
        exit("Error, must specify --genome_lib_dir or set env var CTAT_GENOME_LIB")



    ref_genome_fa = os.path.join(genome_lib_dir, "ref_genome.fa")

    ref_genome_seq_lengths = get_seq_lengths(ref_genome_fa)
    ref_chromosomes = list(filter(lambda x: re.match("chr", x), ref_genome_seq_lengths.keys()))
    

    viral_seq_lengths = get_seq_lengths(virus_db)


    # open some output files
    sim_insertion_targets_ofh = open("{}.insertion_seqs.fa".format(out_prefix), "wt")
    sim_left_fq_ofh = open("{}.insertion_seqs.left_fq".format(out_prefix), "wt")
    sim_right_fq_ofh = open("{}.insertion_seqs.right_fq".format(out_prefix), "wt")

    for virus_acc, seq_len  in viral_seq_lengths.items():
        #print("\t".join([virus_acc, str(seq_len)]))

        for insertion_number in range(1, ins_per_virus+1):


            found = False

            num_tries = 0

            while not found:

                num_tries += 1

                if num_tries > 100:
                    # give up on this one.
                    logger.info("-too many tries for {}, skipping it".format(virus_acc))
                    break
                
                # get random genome seq
                rand_chrom, rand_chrom_lend, rand_chrom_rend, rand_chrom_orient, rand_genome_seq = get_random_seq(ref_genome_fa, ref_chromosomes, ref_genome_seq_lengths, seq_extraction_len)

                # get random viral seq
                rand_virus, rand_virus_lend, rand_virus_rend, rand_virus_orient, rand_virus_seq = get_random_seq(virus_db, [virus_acc], viral_seq_lengths, seq_extraction_len)


                insertion_accession = "~".join([rand_chrom, str(rand_chrom_lend), str(rand_chrom_rend), rand_chrom_orient,
                                                rand_virus, str(rand_virus_lend), str(rand_virus_rend), rand_virus_orient])


                if extreme_N(rand_genome_seq) or extreme_N(rand_virus_seq):
                    logger.info("-skipping {} due to extreme N content".format(insertion_accession))
                    continue

                found = True
                insertion_sequence = rand_genome_seq.upper() + rand_virus_seq.lower()

                insertion_fa_record = ">{}\n{}\n".format(insertion_accession, insertion_sequence)

                logger.info("-processing {}".format(insertion_accession))

                print(insertion_fa_record, file=sim_insertion_targets_ofh, end='')

                # sim fastqs
                simulate_reads(insertion_accession, insertion_sequence, depth_of_coverage, sim_left_fq_ofh, sim_right_fq_ofh, read_length)


    logger.info("done.")
    
    sys.exit(0)



def get_seq_lengths(fasta_filename):

    fai_filename = fasta_filename + ".fai"
    if not os.path.exists(fai_filename):
        cmd = "samtools faidx " + fasta_filename
        subprocess.check_call(cmd, shell=True)

    seq_lengths = dict()
    with open(fai_filename) as fh:
        for line in fh:
            line = line.rstrip()
            vals = line.split("\t")
            acc, length = vals[0], int(vals[1])
            seq_lengths[acc] = length

    return seq_lengths



def extract_sequence(fasta_filename, accession, lend, rend, orient):

    cmd = "samtools faidx {} \'{}\':{}-{}".format(fasta_filename, accession, lend, rend)

    if orient == "-":
        cmd += " --reverse-complement"

    sequence = "".join(subprocess.check_output(cmd, shell=True).decode().rstrip().split("\n")[1:])

    return sequence





def get_random_seq(fasta_filename, accession_list, fasta_seq_lengths, seq_extraction_len):

    random_chrom = random.choice(accession_list)
    chrom_len = fasta_seq_lengths[random_chrom]
    
    selected_chrom_lend = random.randint(1, max(1,chrom_len - seq_extraction_len))
    selected_chrom_rend = min(chrom_len, selected_chrom_lend + seq_extraction_len - 1)
    
    random_seq_orient = random.choice(['+', '-'])

    extracted_seq = extract_sequence(fasta_filename, random_chrom,
                                     selected_chrom_lend, selected_chrom_rend,
                                     random_seq_orient)

    return random_chrom, selected_chrom_lend, selected_chrom_rend, random_seq_orient, extracted_seq


def simulate_reads(accession, sequence, depth, left_fq_ofh, right_fq_ofh, read_length):


    # how many reads do we need to achieve target depth of coverage?
    num_reads = math.ceil(len(sequence) * depth / 2.0 / read_length) 


    tmp_fa = "tmp.fa"
    with open(tmp_fa, "wt") as ofh:
        print(">{}\n{}".format(accession, sequence), file=ofh)

    tmp_left_fq = "tmp.left_fq"
    tmp_right_fq = "tmp.right_fq"
    
    cmd = "wgsim-trans {} -N {} {} {}".format(tmp_fa, num_reads, tmp_left_fq, tmp_right_fq)

    subprocess.check_call(cmd, shell=True)

    # add to our collection of simulated reads.
    with open(tmp_left_fq) as fh:
        for line in fh:
            print(line, file=left_fq_ofh, end='')

    with open(tmp_right_fq) as fh:
        for line in fh:
            print(line, file=right_fq_ofh, end='')


    # cleanup
    os.remove(tmp_fa)
    os.remove(tmp_left_fq)
    os.remove(tmp_right_fq)


    return


def extreme_N(sequence):

    N_count = 0

    seqlen = len(sequence)
    for char in sequence:
        if char.upper() == "N":
            N_count += 1

    if N_count / seqlen > 0.1:
        return True
    else:
        return False

 
####################
 
if __name__ == "__main__":
    main()
