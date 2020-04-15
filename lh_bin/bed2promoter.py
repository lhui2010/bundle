#!/usr/bin/env python
import re
import argparse
from Bio import SeqIO
#import sys

#Args
parser = argparse.ArgumentParser(description='transform a gene bed to promoter bed')
parser.add_argument('INPUT', type=str, nargs = 1,
                    help='input bed file')
parser.add_argument('--version', action='version', version='%(prog)s 0.1')
args = parser.parse_args()

offset = 2000

with open(args.INPUT[0]) as fh:
    for line in fh:
#8_RaGOO	164705426	164705772	A188G13745-t1	0	+
        (chromosome, gene_start, gene_end, gene, score, strand) = line.rstrip().split('\t')
        if(strand == '-'):
            #print('strand - ')
            promoter_start = str( int(gene_end)  + 1 )
            promoter_end = str( int(gene_end) + offset)
        else:
            promoter_start = str(int(gene_start) - offset) if int(gene_start) > offset else '1'
            promoter_end = str( int(gene_start) - 1)


        print ("\t".join([chromosome, promoter_start, promoter_end, gene, score, strand]))

