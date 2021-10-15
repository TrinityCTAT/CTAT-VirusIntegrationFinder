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
import os, re
import logging
import csv
import sys, time
import argparse
import multiprocessing as mp

logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s : %(levelname)s : %(message)s',
                    datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)

#!/usr/bin/env python

import argparse
import math
import string
import os, sys
import logging 
logging.basicConfig(format='\n %(levelname)s : %(message)s', level=logging.DEBUG)
# logging.basicConfig(format='%(asctime)s %(message)s') 


def CHECK_reading(true_ids, read_names):
    if len(true_ids) == len(read_names):
        print(f"\t\t\t\tAll {len(read_names)} reads found!")
    else:
        m = f"\t\t\t\t NOT All reads found! {true_ids2} out of {len(read_names)} Found."
        logger.warning(f"\t\t\t\tAll {len(read_names)} reads found!")
        exit()




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
        # super(faFile, self).__init__()
        # super().__init__()
        self.input_file = input_file
        self.buffer_size = os.stat(input_file).st_blksize
        self.file = open(input_file, "r")
        self.sequence_ids = []
        self.sequences = {}
        self.num_seqs = 0
        self.output_name = ""
        # self.fa_reader()

    def faReader(self):
        contig_id = None
        # Read in the file 
        tmp = self.file.read()
        # Split the file by the > character
        ## this seperates the different sequences
        ## skip the first as it will be blank do to .split() 
        tmp = tmp.split(">")[1:]
        # Create the dictionary to hold the ids and sequences 
        dic = {}
        # Run through each sequence and remove the \n characters 
        ## Then add it to the dictionary 
        for i in tmp:
            idx = i.find("\n")
            dic[i[0:idx]] = i[idx:].replace("\n","")
        # Save the dictionary in the object as .sequences 
        self.sequences = dic
        # return the object
        return self 

    def fqReader(self):
        '''
        Read in the FastQ formated file 
        This function assumes the file is in standard 4 line FASTQ format 

        FastQ file format
            line1: @sequence_identifier 
            line2: Raw_Sequence
            line3: option
            line4: Quality values s
        '''

        contig_id = None
        # Read in the file 
        tmp = self.file.read().rstrip()
        # split all the lines by \n 
        tmp = tmp.split("\n")
        # Create the dictionary to hold the ids and sequences 
        dic = {}
        # Run through each sequence and remove the \n characters 
        ## Then add it to the dictionary 
        n = 4 # there are 4 lines in a single entry
        # for i in list(range(0,len(tmp),n)):
        #     dic[tmp[i][1:]] = "\n".join(tmp[i+1:i+n])+"\n"
        dic = {tmp[i][1:] : "\n".join(tmp[i+1:i+n])+"\n" for i in list(range(0,len(tmp),n))}
        
        # # Save the dictionary in the object as .sequences 
        # self.sequences = dic

        # Convert into Pandas DataFrame for easy subsetting 
        # Need to adjust the ID values as there is a tag on the end 
        df = pd.DataFrame.from_dict(dic, orient = "index")
        df = df.reset_index()
        df.columns = ["ID","Sequence"]
        a = [i.split(" ")[0] for i in df.ID]
        df["ID2"] = a
        
        self.df = df
        
        # return the object
        return self 



    # def faWriter(self, output = "Output_fa.fa"):
    #     # Set the constant 
    #     self.output_name = output

    #     # Create the file to write to
    #     output_file = open(self.output_name, "w")

    #     for i in self.sequences:
    #         tmp = ">{}\n{}\n".format(i, self.sequences[i])
    #         output_file.write(tmp)
    #     output_file.close()



# Create class object faqFile
class ExtractEvidenceReads:
    '''
    Class to extract the reads of interest from given FASTQ file(s) 
    '''
    def __init__(self, args): # arguments to class instantiation 
        
        # Class ExtractEvidenceReads inherits from class faFile
        # super(ExtractEvidenceReads, self).__init__()
        # super().__init__()
        #~~~~~~~~~~~~~~~~~~~~~~~~~
        # Add constants to object 
        #~~~~~~~~~~~~~~~~~~~~~~~~~
        # FASTQs
        self.left_fq  = args.left_fq
        if args.right_fq:
            self.right_fq = args.right_fq
        else: 
            self.right_fq = None
        # Insertion Candidate input 
        self.insertion_candidates = args.insertion_candidates
        # Prefix to use 
        self.out_prefix = args.out_prefix
        


    def parsInsertionCandidates(self):
        '''
        Parse the insertion candidates file 
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
        
        pysam.FastxFile reads through the fastQ files 
        each entry is of class pysam.FastqProxy
        entry == pysam.FastqProxy class
        
        '''
        logger.info("\n################################\n Reading in FASTQ \n################################")
        
        # for entry in fh:
        #     print(entry.name)
        #     print(entry.sequence)
        #     print(entry.comment)
        #     print(entry.quality)

        
        
        total_count = len(self.read_names)
        
        #~~~~~~~~~~~~~~~~
        # Left FastQ file
        #~~~~~~~~~~~~~~~~
        logger.info("\tReading in LEFT FASTQ")

        # Create the FASTQ object 
        fastx_obj = faFile(self.left_fq)
        # Read the FASTQ file into the object 
        fastx_obj = fastx_obj.fqReader()

        # identify which reads to subset 
        true_ids = list(set(self.read_names) & set(fastx_obj.df["ID2"]))
        
        ########
        # CHECK 
        ########
        # Make sure all reads are found 
        CHECK_reading(true_ids,self.read_names)

        # Subset the DF to the reads of interest 
        subset_df = fastx_obj.df[fastx_obj.df['ID2'].isin(true_ids)]
        subset_df = subset_df[["ID","Sequence"]]

        self.df_left = subset_df


        #~~~~~~~~~~~~~~~~
        # Right FastQ file
        #~~~~~~~~~~~~~~~~
        logger.info("\tReading in RIGHT FASTQ")

        if self.right_fq != None:
            # Create the FASTQ object 
            fastx_obj = faFile(self.right_fq)
            # Read the FASTQ file into the object 
            fastx_obj = fastx_obj.fqReader()

            # identify which reads to subset 
            true_ids = list(set(self.read_names) & set(fastx_obj.df["ID2"]))
            
            ########
            # CHECK 
            ########
            # Make sure all reads are found 
            CHECK_reading(true_ids,self.read_names)

            # Subset the DF to the reads of interest 
            subset_df = fastx_obj.df[fastx_obj.df['ID2'].isin(true_ids)]
            subset_df = subset_df[["ID","Sequence"]]

            self.df_right = subset_df

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
        # Create new fastQ file 
        temp_list = []
        for index, row in self.df_left.iterrows():
            # temp_list.append("@" + row["ID"] + "\n" + row["Sequence"])
            temp_list.append( f"@{row['ID']}\n{row['Sequence']}" )
        
        final_str = "".join(temp_list)
        #final_str = final_str.rstrip("\n")
        output_file = open(out_filename1, "w")
        output_file.write(final_str)
        output_file.close()


        #~~~~~~~~~~~~~~~~
        # Right FastQ file
        #~~~~~~~~~~~~~~~~
        if self.right_fq != None:
            logger.info(f"\tWriting Right FASTQ to: {out_filename2}")
            # Create new fastQ file 
            temp_list = []
            for index, row in self.df_right.iterrows():
                # temp_list.append("@" + row["ID"] + "\n" + row["Sequence"])
                temp_list.append( f"@{row['ID']}\n{row['Sequence']}" )
            
            final_str = "".join(temp_list)
            #final_str = final_str.rstrip("\n")
            output_file = open(out_filename2, "w")
            output_file.write(final_str)
            output_file.close()



##########################################
# MAIN
def main():

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Inputs
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    parser = argparse.ArgumentParser(   formatter_class = argparse.RawTextHelpFormatter, 
                                        description     = "")
    parser.add_argument('--left_fq',  required = True,  help = "fastq File input.")
    parser.add_argument('--right_fq', required = False, help = "fastq File input.")
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








