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


class humanChimericAlignments:
    '''
    Class to Define and remove the human-human chimeric alignments from the evidence. 
    '''
    def __init__(self, args): # arguments to class instantiation 

        logger.info("\n################################\n Reading inputs \n################################")
        #~~~~~~~~~~~~~~~~~~~~~~~~~
        # Add constants to object 
        #~~~~~~~~~~~~~~~~~~~~~~~~~
        # Read in the human chimeric reads file 
        logger.info("\tHuman Chimeric Alignment File")
        self.chimeric_df = pd.read_csv(args.human_chimeric_alignments, sep = "\t")
        # Insertion Candidate input 
        logger.info("\tInsertion Candidates")
        self.insertion_candidates = pd.read_csv(args.insertion_candidates, sep = "\t")
        # Prefix to use 
        self.out_prefix = args.out_prefix

    def reviseInsertions(self):
        '''
        Function to revise the current Insertion Candidates. 
            Identifies reads found as human-human and human-virus chimeric reads 
            adjusts the total evidence read count and list of insertion reads for these identified reads.
            Adds an extra column to the data frame holding the ratio of reads removed to total readsa for
            the given insertion. 
        '''

        logger.info("\n################################\n Updating Reads and Read Counts \n################################")

        # Lists that will hold the new columns 
        new_total_column = []
        new_read_support = []
        ratio = []

        # Make a set of the human-human chimeric reads
        human_reads = set(self.chimeric_df["read_name"])

        # Iterate over all insertions 
        for index, row in self.insertion_candidates.iterrows():
            reads = set(row["readnames"].split(","))
            # get the reads in common for the given insertion 
            reads_in_common = list(human_reads & set(reads))

            new_total_column.append(row["total"] - len(reads_in_common))

            # Create the ratio of reads removed and total reads for the given insertion 
            # ratio.append( (row["total"] - len(reads_in_common)) / row["total"])
            ratio.append( (len(reads_in_common)) / row["total"])

            # Remove the human-human chimeric reads form the list of reads 
            new_read_support.append(list(set(reads) - human_reads))

        # Add the new columns to the data frame 
        self.insertion_candidates["readnames"] = new_read_support
        self.insertion_candidates["total"] = new_total_column
        self.insertion_candidates["ratio"] = ratio

        return self


    def outputFile(self):
        logger.info("\n################################\n Writing New Insertion File \n################################")

        output_file = f"{self.out_prefix}.revised.insertion_candidates.tsv"
        self.insertion_candidates.to_csv(output_file, index = False)




##########################################
# MAIN
def main():

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Inputs
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    parser = argparse.ArgumentParser(   formatter_class = argparse.RawTextHelpFormatter, 
                                        description     = "")
    parser.add_argument('--human_chimeric_alignments',  required = True,  help = "STAR Chimeric File input.")
    parser.add_argument('--insertion_candidates', required = True, help = "Names of the evidence reads from which you would extract from the fastq files.")
    parser.add_argument('--out_prefix', required = False, default = "human_chim_pruned", help = "fastq File input ")
    args = parser.parse_args()

    obj = humanChimericAlignments(args)
    print("obj.chimeric_df:\n",obj.chimeric_df)
    print("obj.insertion_candidates:\n",obj.insertion_candidates)

    obj.reviseInsertions()

    obj.outputFile()

    # print("obj.insertion_candidates:\n",obj.insertion_candidates)

if __name__ == "__main__":

    main()








