#!/usr/bin/env python

import argparse
import textwrap
import logging

import re
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

#Prepare fasta
def main():
    prog_name = "align each pair using blat"
    usage = "Another python program"

#use of textwrap allowed custome format like "\t\tusage"
    parser = argparse.ArgumentParser(
        prog=prog_name, 
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(usage), 
        epilog="")
    parser.add_argument("pair", help="pair file of contig IDs(like chr1_A188 chr1_B73")
    parser.add_argument("fasta", help="input fasta(>chr1_A188\nAAACCCT..")
#    parser.add_argument("-f", "--flanking", default=10000, type=int, help="flanking distance default (1000)")
    args = parser.parse_args()  
    pair_file = args.pair
    fasta_file = args.fasta
#Read assembly file
#Read fasta file
    logging.warning("Reading fasta sequences..\n")
    fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    with open(pair_file) as fh:
        tag=False
        for line in fh:
            (qry, ref) = line.rstrip().split()
            with open("qryb", 'w') as out_qry, open("refb", 'w') as out_ref:
                try:
                    out_qry.write(fasta_dict[qry].format('fasta'))
                    out_ref.write(fasta_dict[ref].format('fasta'))
                except KeyError:
                    continue

            #os.system("lastz --gfextend --chain --gapped --ambiguous=n --format=mapping --identity=97  --continuity=1 {} {} > {} ".format('qry', 'ref', "qry.lastz"))
            os.system("blat {} {}  {} ".format('qryb', 'refb', "qry.psl"))
            if(tag):
                os.system("tail -n +6 qry.psl >> total.psl")
            else:
                os.system("tail -n +3 qry.psl >> total.psl")
            tag = True

if __name__ == "__main__":
    main()
