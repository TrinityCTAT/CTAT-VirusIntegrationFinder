#!/usr/bin/env python

import sys, os, re
import pysam
import textwrap



def main():
    usage = "\n\n\tusage: target.fasta kmer_len\n\n"

    if len(sys.argv) < 3:
        exit(usage)


    target_fasta_filename = sys.argv[1]
    kmer_len = int(sys.argv[2])

    all_seen_kmers = set()

    with pysam.FastxFile(target_fasta_filename) as fh:
        for entry in fh:
            #print(entry.name)
            #print(entry.sequence)
            sequence = entry.sequence
            seen_kmer_pos = evaluate_kmers(sequence, kmer_len, all_seen_kmers)

            if seen_kmer_pos:
                sequence = mask_kmers(sequence, seen_kmer_pos, kmer_len)

            print(">{}\n{}".format(entry.name, "\n".join(textwrap.wrap(sequence, 60)).rstrip()))


    sys.exit(0)


def evaluate_kmers(sequence, kmer_len, all_seen_kmers):

    already_seen_kmers = list()

    for i in range(len(sequence)-kmer_len):
        kmer = sequence[i:i+kmer_len]
        #print(kmer)
        kmer_hashcode = hash(kmer)
        if kmer_hashcode in all_seen_kmers:
            already_seen_kmers.append(i)
        all_seen_kmers.add(kmer_hashcode)

    return already_seen_kmers


def mask_kmers(sequence, kmer_pos_list, kmer_len):

    last_masked_pos = -1

    sequence = list(sequence)
    
    for kmer_start in kmer_pos_list:
        mask_begin_pos = kmer_start
        mask_end_pos = kmer_start + kmer_len
        if last_masked_pos > 0:
            mask_begin_pos = max(mask_begin_pos, last_masked_pos + 1)
            
        for pos in range(mask_begin_pos, mask_end_pos):
            sequence[pos] = 'N'

        last_masked_pos = mask_end_pos - 1


    sequence = "".join(sequence)

    return sequence



if __name__=='__main__':
    main()
