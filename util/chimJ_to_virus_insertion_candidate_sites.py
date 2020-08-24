#!/usr/bin/env python3

import sys, os, re
import argparse
import subprocess
import logging
from collections import defaultdict

if sys.version_info[0] < 3:
    raise Exception("This script requires Python 3")

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s : %(levelname)s : %(message)s',
                    datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)


'''  Chimeric.out.junction file formatting from the STAR manual:

        # from star doc:
        #The rst 9 columns give information about the chimeric junction:

        #The format of this le is as follows. Every line contains one chimerically aligned read, e.g.:
        #chr22 23632601 + chr9 133729450 + 1 0 0 SINATRA-0006:3:3:6387:56650 23632554 47M29S 133729451 47S29M40p76M
        #The first 9 columns give information about the chimeric junction:

        #column 1: chromosome of the donor
        #column 2: rst base of the intron of the donor (1-based)
        #column 3: strand of the donor
        #column 4: chromosome of the acceptor
        #column 5: rst base of the intron of the acceptor (1-based)
        #column 6: strand of the acceptor
        #column 7: junction type: -1=encompassing junction (between the mates), 1=GT/AG, 2=CT/AC
        #column 8: repeat length to the left of the junction
        #column 9: repeat length to the right of the junction
        #Columns 10-14 describe the alignments of the two chimeric segments, it is SAM like. Alignments are given with respect to the (+)
 strand
        #column 10: read name
        #column 11: rst base of the rst segment (on the + strand)
        #column 12: CIGAR of the rst segment
        #column 13: rst base of the second segment
        #column 14: CIGAR of the second segment
'''



def main():

    parser = argparse.ArgumentParser(description="Defines candidate virus insertion sites based on chimeric reads",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("--chimJ", type=str, required=True, help="STAR Chimeric.out.junction file")
    
    parser.add_argument("--patch_db_fasta", type=str, required=True,
                        help="database of the additional genome targets")


    args_parsed = parser.parse_args()

    chimJ_filename = args_parsed.chimJ
    patch_db_fasta_filename = args_parsed.patch_db_fasta


    ## get list of patch_db entries.
    patch_db_entries = set()
    with open(patch_db_fasta_filename, 'rt') as fh:
        for line in fh:
            m = re.search("^>(\S+)", line)
            if m:
                acc = m.group(1)
                patch_db_entries.add(acc)

    
    junction_type_encoding = { "-1" : "Span",
                              "0" : "Split",
                              "1" : "GT/AG",
                              "2" : "CT/AC" }

    opposite_orientation = { '+' : '-',
                             '-' : '+' }


    genome_pair_to_evidence = defaultdict(set)


    duplicate_read_catcher = set()
    
    with open(chimJ_filename, 'rt') as fh:
        for line in fh:
            line = line.rstrip()
            vals = line.split("\t")
            junction_type = junction_type_encoding[ vals[6] ]
            readname = vals[9]
            (chrA, coordA, orientA) = (vals[0], vals[1], vals[2]);
            (chrB, coordB, orientB) = (vals[3], vals[4], vals[5]);

            if not (chrA in patch_db_entries) ^ (chrB in patch_db_entries):
                # must involve just the host genome and a patch db entry
                continue
            
            ## reorient so host genome is always in the + reference orientation.
            if (
                (chrA in patch_db_entries and  orientB == '-')
                or
                (chrB in patch_db_entries and orientA == '-') ):


                token = hash("^".join(vals[0:6] + vals[10:]))
                if token in duplicate_read_catcher:
                    continue

                duplicate_read_catcher.add(token)
                
                # swap them so that the host chr shows up in + orient.
                (chrA, coordA, orientA,
                 chrB, coordB, orientB) = (chrB, coordB, opposite_orientation[orientB],
                                           chrA, coordA, opposite_orientation[orientA] )

            
            chim_read = Chimeric_read(chrA, coordA, orientA, chrB, coordB, orientB, junction_type)

            genome_pair = "^".join([chrA, chrB])

            genome_pair_to_evidence[genome_pair].add(chim_read)


    for genome_pair in genome_pair_to_evidence:
        print(genome_pair)
        for chim_read in genome_pair_to_evidence[genome_pair]:
            print("\t" + str(chim_read))

    
    sys.exit(0)




class Chimeric_read:

    def __init__(self, chrA, coordA, orientA, chrB, coordB, orientB, splitType):
        self.chrA = chrA
        self.coordA = coordA
        self.orientA = orientA
        self.chrB = chrB
        self.coordB = coordB
        self.orientB = orientB
        self.splitType = splitType

    def __repr__(self):
        return "\t".join([self.chrA, self.coordA, self.orientA, self.chrB, self.coordB, self.orientB, self.splitType])

    






if __name__ == '__main__':
    main()
