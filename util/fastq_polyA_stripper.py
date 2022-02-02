#!/usr/bin/env python3
# encoding: utf-8

import os, sys, re
import logging
import argparse
import gzip
from collections import defaultdict

def main():

    parser = argparse.ArgumentParser(description="strips polyA from reads", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--left_fq", required=True, type=str, help="left fastq file")

    parser.add_argument("--right_fq", required=False, type=str, help="right fastq file")

    parser.add_argument("--min_strip_len", default = 10, help='min length to strip from terminal of sequence with poly-A')

    parser.add_argument("--match_score", default=1, help='match score')
    parser.add_argument("--mismatch_penalty", default=-2, help='mismatch penalty')
    parser.add_argument("--max_drop", default=-8, help='max drop in score to cease search')
    parser.add_argument("--min_seq_len", default=30, help="minimum sequence length after trimming")

    
    args = parser.parse_args()

    left_fq = args.left_fq
    right_fq = args.right_fq
    min_strip_len = args.min_strip_len
    match_score = args.match_score
    mismatch_penalty = args.mismatch_penalty
    max_drop = args.max_drop
    min_seq_len = args.min_seq_len

    trim_config = { 'min_strip_len' : min_strip_len,
                    'match_score' : match_score,
                    'mismatch_penalty' : mismatch_penalty,
                    'max_drop' : max_drop,
                    'min_seq_len' : min_seq_len
                    }
    

    left_fq_iterator = fastq_iterator(left_fq)
    right_fq_iterator = None
    if right_fq:
        right_fq_iterator = fastq_iterator(right_fq)


    left_trimmed_outfile = os.path.basename(left_fq) + ".polyA-trimmed.fastq"
    left_ofh = open(left_trimmed_outfile, "wt")
    
    if right_fq:
        right_trimmed_outfile = os.path.basename(right_fq) + ".polyA-trimmed.fastq"
        right_ofh = open(right_trimmed_outfile, "wt")
    

    counter = defaultdict(int)
    
    for left_read_tuple in left_fq_iterator:
        left_readname, left_readseq, left_L3, left_quals = left_read_tuple
        left_readseq_len = len(left_readseq)
        
        left_readseq, left_quals = polyA_trim(left_readseq, left_quals, trim_config)


        if left_readseq != "":
            if right_fq_iterator is not None:
                right_read_tuple = next(right_fq_iterator)
                right_readname, right_readseq, right_L3, right_quals = right_read_tuple
                right_readseq_len = len(right_readseq)

                right_readseq, right_quals = polyA_trim(right_readseq, right_quals, trim_config)
                # require both pairs to pass trimming requirements
                if right_readseq != "":
                    print("\n".join([left_readname, left_readseq, left_L3, left_quals]), file=left_ofh)
                    print("\n".join([right_readname, right_readseq, right_L3, right_quals]), file=right_ofh)

                    if len(left_readseq) != left_readseq_len or len(right_readseq) != right_readseq_len:
                        counter['PE_TRIMmed_reported'] += 1
                    else:
                        counter['PE_untrimmed_reported'] += 1
                        
                else:
                    counter['PE_rejected'] += 1

            else:
                # unpaired, just process left fq
                print("\n".join([left_readname, left_readseq, left_L3, left_quals]), file=left_ofh)
                if len(left_readseq) != left_readseq_len:
                    counter['SE_TRIMmed_reported'] += 1
                else:
                    counter['SE_untrimmed_reported'] += 1

        else:
            counter['SE_rejected'] += 1


    print(str(counter))
    
    sys.exit(0)




def polyA_trim(readseq, quals, trim_config):

    # trim terminal polyA
    readseq, quals = terminal_trim(readseq, quals, 'A', trim_config)

    if readseq:
        readseq, quals = terminal_trim(readseq, quals, readseq[-1], trim_config)
    
        
    # trim initial polyT
    if readseq:
        readseq, quals = initial_trim(readseq, quals, 'T', trim_config)

    if readseq:
        readseq, quals = initial_trim(readseq, quals, readseq[0], trim_config)
    
    if len(readseq) >= trim_config['min_seq_len']:
        return(readseq, quals)
    else:
        return("", "")




def terminal_trim(readseq, quals, nucleotide, trim_config):
    trim_rend_pos = compute_trim_pos(readseq, nucleotide, trim_config)
    if trim_rend_pos is not None:
        readseq = readseq[:trim_rend_pos]
        quals = quals[:trim_rend_pos]
        
    return readseq, quals



def initial_trim(readseq, quals, nucleotide, trim_config):
    
    revreadseq = readseq[::-1]
    trim_pos = compute_trim_pos(revreadseq, nucleotide, trim_config)
    if trim_pos is not None:
        # revcomp the trim pos
        trim_pos = len(revreadseq) - trim_pos - 1
        readseq = readseq[trim_pos+1:]
        quals = quals[trim_pos+1:]
        
    return readseq, quals


def compute_trim_pos(readseq, nucleotide, trim_config):

    score = 0
    max_score = 0
    max_score_pos = -1
    
    readseq = readseq.upper()
    for i in range(len(readseq)-1, -1, -1):
        if readseq[i] == nucleotide:
            score += trim_config['match_score']
        else:
            score += trim_config['mismatch_penalty']

        if score > max_score:
            max_score = score
            max_score_pos = i

        if score - max_score <= trim_config['max_drop']:
            break

    if max_score >= trim_config['min_strip_len'] * trim_config['match_score']:
        return max_score_pos
    else:
        return None
    

def fastq_iterator(fastq_filename):

    if re.search(".gz$", fastq_filename):
        fh = gzip.open(fastq_filename, 'rt')
    else:
        fh = open(fastq_filename, 'rt')

    have_records = True
    while have_records:
        readname = next(fh).rstrip()
        readseq = next(fh).rstrip()
        L3 = next(fh).rstrip()
        quals = next(fh).rstrip()

        yield (readname, readseq, L3, quals)

        if not readname:
            break

    return
        



####### unit testing

def test_trim_poly():

    trim_config = { 'min_strip_len' : 10,
                    'match_score' : 1,
                    'mismatch_penalty' : -2,
                    'max_drop' : -6,
                    'min_seq_len' : 5
                    }

    readseq = "TTTTTTTTTTTTTTTTGATCGATCGATCAAAAAAAAAAAAAAA"
    trimmed_seq, trimmed_quals = polyA_trim(readseq, "2" * len(readseq), trim_config)

    assert trimmed_seq == "GATCGATCGATC", f"Error {trimmed_seq} not as expected" 




if __name__=='__main__':
    main()
