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


    parser.add_argument("--aggregation_dist", type=int, required=False, default=100,
                        help="distance around top chimeric event breakpoint for aggregating supporting reads")

    parser.add_argument("--output_prefix", "-o", dest="output_prefix", type=str, required=True,
                        help = "output prefix")
    
    args_parsed = parser.parse_args()

    chimJ_filename = args_parsed.chimJ
    patch_db_fasta_filename = args_parsed.patch_db_fasta
    aggregation_dist = args_parsed.aggregation_dist
    output_prefix = args_parsed.output_prefix
    
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
                               "1" : "Split", #GT/AG",
                               "2" : "Split" } #CT/AC" }

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

            coordA = int(coordA)
            coordB = int(coordB)
            
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

            
            chim_read = Chimeric_read(chrA, coordA, orientA, chrB, coordB, orientB, junction_type, readname)

            genome_pair = "^".join([chrA, chrB])

            genome_pair_to_evidence[genome_pair].add(chim_read)


    all_chim_events = list()

    for genome_pair in genome_pair_to_evidence:

        chim_reads = genome_pair_to_evidence[genome_pair]

        logger.info("Genome pair: " + genome_pair)

        chim_events = group_chim_reads_into_events(chim_reads, aggregation_dist)

        # //TODO: can try to link up events here.  

        all_chim_events.extend(chim_events)


    ## prioritize by total read support.
    all_chim_events = sorted(all_chim_events, key=lambda x: x.get_read_support()[2], reverse=True)
    

    ## generate final report.
    output_filename_full = output_prefix + ".full.tsv"
    output_filename_abridged = output_prefix + ".abridged.tsv"
    with open(output_filename_abridged, 'wt') as ofh:
        with open(output_filename_full, 'wt') as ofh_full:

            # print header
            ofh.write("\t".join(["entry", "chrA", "coordA", "orientA", "chrB", "coordB", "orientB", "primary_brkpt_type", "num_primary_reads", "num_supp_reads", "total_reads"]) + "\n") 
            ofh_full.write("\t".join(["entry", "chrA", "coordA", "orientA", "chrB", "coordB", "orientB", "primary_brkpt_type", "num_primary_reads", "num_supp_reads", "total_reads", "readnames"]) + "\n") 
            
            chim_counter = 0
            for chim_event in all_chim_events:
                chim_counter += 1
                print(str(chim_counter) + "\t" + str(chim_event), file=ofh)
                supporting_reads = chim_event.get_readnames()
                print(str(chim_counter) + "\t" + str(chim_event) + "\t" + ",".join(supporting_reads), file=ofh_full)

    logger.info("-wrote output to {}".format(output_filename_abridged))
    
    sys.exit(0)




def group_chim_reads_into_events(chim_reads_list, aggregation_dist):

    chim_events = list()

    remaining_reads = chim_reads_list

    while remaining_reads:

        top_event_reads, remaining_reads = gather_top_event_reads(remaining_reads)
        chim_event = Chimeric_event(top_event_reads)
        
        if not supplements_existing_event(chim_event, chim_events, aggregation_dist):
            chim_events.append(chim_event)
            logger.info('-logging chimeric event: ' + str(chim_event))

    return chim_events

    
def supplements_existing_event(chim_event, chim_events_list, aggregation_dist):

    for prev_chim_event in chim_events_list:

        if ( chim_event.orientA == prev_chim_event.orientA
             and
             abs(chim_event.coordA - prev_chim_event.coordA) <= aggregation_dist
             and
             chim_event.orientB == prev_chim_event.orientB
             and
             abs(chim_event.coordB - prev_chim_event.coordB) <= aggregation_dist ):

            logger.info("-adding {} as supplement to {}".format(str(chim_event), str(prev_chim_event)))

            prev_chim_event.absorb_supporting_reads(chim_event.chimeric_reads_list)

            return True

    return False


        



def gather_top_event_reads(reads_list):

    # count reads according to breakpoint.

    brkpt_counter = defaultdict(int)
    brkpt_type = dict()
    
    for read in reads_list:
        brkpt = "{}^{}".format(read.coordA, read.coordB)
        brkpt_counter[brkpt] += 1
        brkpt_type[brkpt] = read.splitType
    
    # prioritize split reads over spanning reads.
    priority = { 'Split' : 1,
                 'Span' : 0 }

    sorted_brkpts = sorted(brkpt_counter.keys(), key=lambda x: (priority[brkpt_type[x]], brkpt_counter[x]), reverse=True)

    top_brkpt = sorted_brkpts[0]

    top_event_reads = list()
    remaining_reads = list()

    for read in reads_list:
        brkpt = "{}^{}".format(read.coordA, read.coordB)
        if brkpt == top_brkpt:
            top_event_reads.append(read)
        else:
            remaining_reads.append(read)

    return top_event_reads, remaining_reads
    







class Chimeric_read:

    def __init__(self, chrA, coordA, orientA, chrB, coordB, orientB, splitType, readname=None):
        self.chrA = chrA
        self.coordA = coordA
        self.orientA = orientA
        self.chrB = chrB
        self.coordB = coordB
        self.orientB = orientB
        self.splitType = splitType
        self.readname = readname
    
        
    def __repr__(self):
        return "\t".join([self.chrA, str(self.coordA), self.orientA, self.chrB, str(self.coordB), self.orientB, self.splitType])

    

class Chimeric_event (Chimeric_read):


    def __init__(self, chimeric_reads_list):
        self.chimeric_reads_list = chimeric_reads_list

        self.chimeric_reads_absorbed = list()  # for reads that also support this event but are either spanning or split w/ different nearby brkpt

        example_read = chimeric_reads_list[0]
        super().__init__(example_read.chrA, example_read.coordA, example_read.orientA,
                       example_read.chrB, example_read.coordB, example_read.orientB,
                       example_read.splitType)

        
    
    def absorb_supporting_reads(self, chim_reads_list):
        self.chimeric_reads_absorbed.extend(chim_reads_list)


    def get_read_support(self):
        num_chimeric_reads = len(self.chimeric_reads_list)
        num_absorbed_reads = len(self.chimeric_reads_absorbed)
        num_total_reads = num_chimeric_reads + num_absorbed_reads

        return num_chimeric_reads, num_absorbed_reads, num_total_reads

    def __repr__(self):
        num_chimeric_reads, num_absorbed_reads, num_total_reads = self.get_read_support()
        return(super().__repr__() + "\t{}\t{}\t{}".format(num_chimeric_reads, num_absorbed_reads, num_total_reads))
    
    def get_readnames(self):
        readnames = list()
        for chim_read in self.chimeric_reads_list + self.chimeric_reads_absorbed:
            readnames.append(chim_read.readname)

        return readnames



if __name__ == '__main__':
    main()
