#!/usr/bin/env python

import argparse
import textwrap
import re

def detect_telomere(seq):
    motif = [
    re.compile(r'(AACCCTA){5}', re.IGNORECASE), 
    re.compile(r'(ACCCTAA){5}', re.IGNORECASE), 
    re.compile(r'(CCCTAAA){5}', re.IGNORECASE), 
    re.compile(r'(CCTAAAC){5}', re.IGNORECASE),
    re.compile(r'(CTAAACC){5}', re.IGNORECASE),
    re.compile(r'(TAAACCC){5}', re.IGNORECASE),
    re.compile(r'(AAACCCT){5}', re.IGNORECASE),
    re.compile(r'(TAGGGTT){5}', re.IGNORECASE),
    re.compile(r'(AGGGTTT){5}', re.IGNORECASE),
    re.compile(r'(GGGTTTA){5}', re.IGNORECASE),
    re.compile(r'(GGTTTAG){5}', re.IGNORECASE),
    re.compile(r'(GTTTAGG){5}', re.IGNORECASE),
    re.compile(r'(TTTAGGG){5}', re.IGNORECASE),
    re.compile(r'(TTAGGGT){5}', re.IGNORECASE)]

    repeat_seq_len = len('AAACCCT')
    repeat_times_cutoff = 6

    target = ''
    max_len = 0
    repeat_num = 0
    for m in motif:
        res_list = re.search(m, seq)
        if(res_list != None):
            span = re.sub(r'\).*', ')', str(res_list)[29:])
            motif = str(m)[13:20]
            return motif+span
    return 'NA'



def main():
    prog_name = "Another python program"
    usage = "Another python program"

    parser = argparse.ArgumentParser(
        prog=prog_name, 
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(usage), 
        epilog="")
    parser.add_argument("qry1", help="qry1 file")
#    parser.add_argument("-f", "--flanking", default=10000, type=int, help="flanking distance default (1000)")
    args = parser.parse_args()  

    qry1_file = args.qry1
#    flanking_distance = args.flanking
    with open(qry1_file) as fh:
        for line in fh:
            (fasta_name, left_seq, right_seq) = line.rstrip().split()
            left = detect_telomere(left_seq)
            right = detect_telomere(right_seq)
            print("\t".join([fasta_name, left, right]))

if __name__ == "__main__":
    main()
