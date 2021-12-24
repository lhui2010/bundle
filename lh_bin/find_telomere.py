#!/usr/bin/env python

import argparse
import textwrap
import re

def detect_telomere(seq):
    motif = [
    re.compile(r'(AACCCTA){10}', re.IGNORECASE), 
    re.compile(r'(ACCCTAA){10}', re.IGNORECASE), 
    re.compile(r'(CCCTAAA){10}', re.IGNORECASE), 
    re.compile(r'(CCTAAAC){10}', re.IGNORECASE),
    re.compile(r'(CTAAACC){10}', re.IGNORECASE),
    re.compile(r'(TAAACCC){10}', re.IGNORECASE),
    re.compile(r'(AAACCCT){10}', re.IGNORECASE),
    re.compile(r'(TAGGGTT){10}', re.IGNORECASE),
    re.compile(r'(AGGGTTT){10}', re.IGNORECASE),
    re.compile(r'(GGGTTTA){10}', re.IGNORECASE),
    re.compile(r'(GGTTTAG){10}', re.IGNORECASE),
    re.compile(r'(GTTTAGG){10}', re.IGNORECASE),
    re.compile(r'(TTTAGGG){10}', re.IGNORECASE),
    re.compile(r'(TTAGGGT){10}', re.IGNORECASE)]

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
    usage = """fasta_ends.py -f 1000 ${REF} > ${REF}.ends
find_telomere.py ${REF}.ends > ${REF}.ends.telo
grep '(' ${REF}.ends.telo  |wc -l > ${REF}.ends.telo.count
"""

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
