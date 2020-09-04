#!/usr/bin/env python

import argparse
import textwrap
from Bio import SeqIO

def main():
    prog_name = "Another python program"
    usage = "Another python program"

    parser = argparse.ArgumentParser(
        prog=prog_name, 
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(usage), 
        epilog="")
    parser.add_argument("qry1", help="qry1 file")
    parser.add_argument("-f", "--flanking", default=100, type=int, help="flanking distance default (100)")
    args = parser.parse_args()  

    qry1_file = args.qry1
    chunk_size = args.flanking
    rend = chunk_size * -1
    for record in SeqIO.parse(qry1_file, "fasta"):
        print("%s\t%s\t%s" % (record.id, record.seq[:chunk_size], 
            record.seq[rend:]))

if __name__ == "__main__":
    main()
