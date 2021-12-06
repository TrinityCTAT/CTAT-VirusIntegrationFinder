#!/usr/bin/env python3


"""

Go to http://www.virusite.org/

download group=human viruses as 'human_viruses.list.csv'

download all viruses from donwload site, unzip, called 'genomes.fasta'

use this script to extract the corresponding set of human viruses from genomes.fasta:

    extract_human_viruses.py | sed -e 's/,_complete_sequence//' | sed -e 's/,_complete_genome//' >  human_viruses.fasta


"""


import sys, os, re
import pysam
import textwrap



accs_want = set()
with open("human_viruses.list.csv") as fh:
    for line in fh:
        line = line.rstrip()
        if line == "":
            continue
        
        fields = line.split(";")
        acc_val = fields[5]
        print(acc_val)
        accs_want.add(acc_val)


with pysam.FastxFile("genomes.fasta") as fh:
    for entry in fh:
        name = entry.name
        fullname = entry.name + " " + entry.comment
        acc = name.split("|")[1]
        if acc in accs_want:
            sequence = entry.sequence
            sequence = textwrap.wrap(sequence, 60)
            sequence = "\n".join(sequence).rstrip()
            fullname = fullname.replace(" ", "_")

            print(">{}\n{}".format(fullname, sequence))


sys.exit(0)


        
