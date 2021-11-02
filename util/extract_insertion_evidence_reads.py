#!/usr/bin/env python3
# -*- coding: utf-8 -*-


####################
# Preliminary 
####################
#~~~~~~~~~~~~~~~~~~
# import packages 
#~~~~~~~~~~~~~~~~~~
import pandas as pd
import pysam 
import os, sys, re
import logging
import csv
import sys, time
import argparse
import multiprocessing as mp
import gzip
import itertools

logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s : %(levelname)s : %(message)s',
                    datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)



def CHECK_reading(true_ids, read_names):
    if len(true_ids) == len(read_names):
        print(f"\t\t\t\tAll {len(read_names)} reads found!")
    else:
        m = f"\t\t\t\t NOT All reads found! {len(true_ids)} out of {len(read_names)} Found."
        logger.warning(m)

        print(set(read_names) - set(true_ids))
        exit()

def readME(infile, read_names):
    '''
    Read in the fastQ and iterate over it 
    '''
    total_reads = len(read_names)
    final_list = []
    while True:
        i = list(itertools.islice(infile, 4))
        if not i:
            break 
        ID = i[0][1:].split(" ")[0]
        if ID in read_names:
            final_list.append(i)

            # print Percent
            sys.stdout.write("\r"); sys.stdout.flush()
            sys.stdout.write(f"{len(final_list)/total_reads * 100}%...")
    return final_list

def readMeOldFormat(infile, read_names):
    '''
    If given a fastq file with older formating, 
    need to adjust the IDs by removing the "/1" or "/2" at the end 
    '''
    final_list = []
    # Iterate over the fastq file, 4 lines at a time
    # returns the 4 lines as a single list 
    while True:
        i = list(itertools.islice(infile, 4))
        if not i:
            break
        ID = i[0][1:].split(" ")[0][:-2]
        if ID in read_names:
            final_list.append(i)
            
            # print Percent 
            sys.stdout.write("\r"); sys.stdout.flush()
            sys.stdout.write(f"{len(final_list)/total_reads * 100}%...")

    return final_list


# Create class object faqFile
class faFile:
    '''
    Class to read in FastQ files and parse them. 
        Returns the FASTQ file in a pandas DataFrame format for easy parsing 
    
    Input Files:
        FASTA : stores sequence records. stores an ID and a sequence. 
        FASTQ : contains a sequence of quality scores for each nucleotide

    Python: 
        self            : self represents the instance of the class, used to access attributes and methods of the class. maps attributes with given arguments 
        init            : 'init' is a constructor, method called when object created allows for teh initialization of the attributes of a class 
    '''

    # Class ExtractEvidenceReads inherits from class faFile
    def __init__(self, input_file): # arguments to class instantiation 

        # Apply constants to the object 
        self.input_file = input_file
        self.buffer_size = os.stat(input_file).st_blksize
        self.sequence_ids = []
        self.sequences = {}
        self.num_seqs = 0
        self.output_name = ""

        #~~~~~~~~~~~
        # Read File 
        #~~~~~~~~~~~
        # # check if the file is gzipped 
        # if input_file.endswith(".gz"):
        #     self.file = gzip.open(input_file).read().rstrip().decode()
        # else:
        #     self.file = open(input_file, "r").read().rstrip()


    def fqReader(self, read_names):
        '''
        Subsets the given Fastq files to only include the given read_names.
        Returns the new fastq as a sting.

        This function assumes the file is in standard 4 line FASTQ format 
        FastQ file format
            line1: @sequence_identifier 
            line2: Raw_Sequence
            line3: option
            line4: Quality values s
        '''

        # variable to tell if the files are old format or not
        older_format = False

        # If files are gzipped
        if self.input_file.endswith(".gz"):
            with gzip.open(self.input_file, "rt") as infile:
                # check the first line 
                # Need to determine if these are older or newer formatted fastqs
                first_line = infile.readline()
                first_id = first_line[1:].split(" ")[0]
                if ( (first_id.endswith("/1")) or (first_id.endswith("/2")) ):
                    logger.info("\t\tIdentified older formated Fastq, editing read names...")
                    older_format = True
            # Now iterate over file 
            with gzip.open(self.input_file, "rt") as infile:
                if older_format == True:
                    final_output = readMeOldFormat(infile,read_names)
                else:
                    final_output = readME(infile, read_names)

        else:
            with open(self.input_file, "r") as infile:
                # check the first line 
                first_line = infile.readline()
                first_id = first_line[1:].split(" ")[0]
                if ( (first_id.endswith("/1")) or (first_id.endswith("/2")) ):
                    logger.info("\t\tIdentified older formated Fastq, editing read names...")
                    older_format = True

            with open(self.input_file, "r") as infile:
                if older_format == True:
                    final_output = readMeOldFormat(infile,read_names)
                else:
                    final_output = readME(infile, read_names)

        #~~~~~~~~~~~~~~~~~~~~~~~~~
        # CHECK
        #~~~~~~~~~~~~~~~~~~~~~~~~~
        N_reads = len(final_output)
        if N_reads == len(read_names):
            m = f"All Reads Found {N_reads}/{len(read_names)}"
            logger.info(m)
        else:
            m = f"NOT All Reads Found {N_reads}/{len(read_names)}"
            logger.warning(m)
            quit()


        # Now join the list of lists together 
        final_output_join = list(itertools.chain.from_iterable(final_output))

        # now join the string in the list together 
        final_output_joined = "".join(final_output_join)
        

        return final_output_joined



# Create class object faqFile
class ExtractEvidenceReads:
    '''
    Class to extract the reads of interest from given FASTQ file(s) 
    '''
    def __init__(self, args): # arguments to class instantiation 
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~
        # Add constants to object 
        #~~~~~~~~~~~~~~~~~~~~~~~~~
        # FASTQs

        self.fastqs  = args.fastqs

        #~~~~~~~~~~~~~~~~~~~
        # check if directory, the direcortry should hold fastq files 
        #~~~~~~~~~~~~~~~~~~~
        if os.path.isdir(self.fastqs[0]):
            fastq_directory = self.fastqs[0]
            # List files in directory 
            fastq_list = os.listdir(self.fastqs[0])
            self.left_fq = os.path.join(fastq_directory, fastq_list[0])

            logger.info("FASTQ input Files:")
            [print(f"\t\tfastq: {i}") for i in fastq_list]
            # check if there is a second file 
            if len(fastq_list) == 2:
                self.right_fq = os.path.join(fastq_directory, fastq_list[1])
            else:
                self.right_fq = None

            # check if there are more then 2 files in the directory 
            if len(fastq_list) > 2:
                logger.warning("More than two files found in fastq directory!")
                exit()

        # ~~~~~~~~~~~
        # files 
        # ~~~~~~~~~~~
        else:
            # set the left fastq file in the object 
            self.left_fq = self.fastqs[0]
            print(f"\t\tfastq: {self.left_fq}")
            if len(self.fastqs) > 1:
                self.right_fq = self.fastqs[1]
                print(f"\t\tfastq: {self.right_fq}")
            else:
                self.right_fq = None

        # Insertion Candidate input 
        self.insertion_candidates = args.insertion_candidates
        # Prefix to use 
        self.out_prefix = args.out_prefix
        


    def parsInsertionCandidates(self):
        '''
        Parse the insertion candidates file 

        Adds the read names of interest found from the insertion candidates file
            and adds them to the object as a list.
        '''
        # Read in the insertion candidates file 
        logger.info("\n################################\n Reading in Insertion Candidates \n################################")
        df = pd.read_csv(self.insertion_candidates, sep = "\t")
        self.insertion_candidates = df
        # print(df)

        # Get the read names of interest 
        read_name_list = []
        for i in df["readnames"]:
            for j in i.split(","):
                read_name_list.append(j)

        logger.info(f"\tNumber of Reads: {len(read_name_list)}")
        self.read_names = read_name_list

        return self

    def subsetFastqFile(self):
        '''
        Read through the fastQ files, extract the reads of interest and write to a new fastQ file. 
        input
            self : ExtractEvidenceReads object 

        output
            Fastq1 : Left FastQ, subset of the original Left FastQ
            Fastq2 : Right FastQ, subset of the original Right FastQ
        
        '''
        logger.info("\n################################\n Reading in FASTQ \n################################")
        
        
        total_count = len(self.read_names)
        
        #~~~~~~~~~~~~~~~~
        # Left FastQ file
        #~~~~~~~~~~~~~~~~
        logger.info("\tReading in LEFT FASTQ")

        # Create the FASTQ object 
        fastx_obj = faFile(self.left_fq)
        
        # Read the FASTQ file into the object 
        fastq_str = fastx_obj.fqReader(self.read_names)
        

        self.str_left = fastq_str


        #~~~~~~~~~~~~~~~~
        # Right FastQ file
        #~~~~~~~~~~~~~~~~

        if self.right_fq != None:
            logger.info("\tReading in RIGHT FASTQ")
            # Create the FASTQ object 
            fastx_obj = faFile(self.right_fq)
            # Read the FASTQ file into the object 
            fastq_str = fastx_obj.fqReader(self.read_names)

            self.str_right = fastq_str

        return self


    def writeFastq(self):


        logger.info("\n################################\n Writing FASTQ \n################################")
        #~~~~~~~~~~~~~~~
        # Preliminary
        #~~~~~~~~~~~~~~~
        # Specify the output file names 
        out_filename1 = f"{self.out_prefix}_1.fastq"
        out_filename2 = f"{self.out_prefix}_2.fastq"

        #~~~~~~~~~~~~~~~~
        # Left FastQ file
        #~~~~~~~~~~~~~~~~
        logger.info(f"\tWriting Left FASTQ to: {out_filename1}")

        output_file = open(out_filename1, "w")
        output_file.write(self.str_left)
        output_file.close()


        #~~~~~~~~~~~~~~~~
        # Right FastQ file
        #~~~~~~~~~~~~~~~~
        if self.right_fq != None:
            logger.info(f"\tWriting Right FASTQ to: {out_filename2}")
            # Create new fastQ file 
            output_file = open(out_filename2, "w")
            output_file.write(self.str_right)
            output_file.close()



##########################################
# MAIN
def main():

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Inputs
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    parser = argparse.ArgumentParser(   formatter_class = argparse.RawTextHelpFormatter, 
                                        description     = "")
    parser.add_argument('--fastqs',  required = True,  nargs='+', help = "fastq File input, or a directory holding one or two.")
    parser.add_argument('--insertion_candidates', required = True, help = "Names of the evidence reads from which you would extract from the fastq files.")
    parser.add_argument('--out_prefix', required = False, default = "ev_reads", help = "fastq File input ")
    args = parser.parse_args()


    # Create the object 
    obj = ExtractEvidenceReads(args)
    # Get the Read names of interest 
    obj.parsInsertionCandidates()
    # Extract the reads from the FastQ files 
    obj.subsetFastqFile()
    # Now write the new FastQ files 
    obj.writeFastq()


if __name__ == "__main__":

    main()








