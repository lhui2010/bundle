#!/usr/bin/env python
import logging
import os
import sys
import argparse

#TODO
#Remove global variable

#OptParser
usage = """Calculate Identity from orthologs 
Usage: 
    {} Ortho_File diff.table  >diff.table.with_rank
""".format(__file__)

parser = argparse.ArgumentParser()
parser.add_argument("file_sorted_rank", help="sorted rank file created by compute_ortho_rank.py")

args = parser.parse_args()  

sorted_rank_file = args.file_sorted_rank

is_block_transposition = "N"

with open(sorted_rank_file) as fh:
    l = fh.readline()
    mylist = l.split()
    (qry_gene, ref_gene, rank_diff) = mylist[0:3]
    last_rank_diff = 0
    last_l = "\t".join([qry_gene, ref_gene, rank_diff])
    l = fh.readline()
    while(l):
        mylist = l.split()
        (qry_gene, ref_gene, rank_diff) = mylist[0:3]
        if(rank_diff == last_rank_diff):
            is_block_transposition = 'Y'
            print(last_l.rstrip() + "\t" + is_block_transposition)
        elif(is_block_transposition == "Y"):
            print(last_l.rstrip() + "\t" + is_block_transposition)
            is_block_transposition = 'N'
        else:
            is_block_transposition = 'N'
            print(last_l.rstrip() + "\t" + is_block_transposition)
        last_rank_diff = rank_diff
        last_l = "\t".join(mylist[0:3])
        l = fh.readline()
        if(l == ""):
            print(last_l.rstrip() + "\t" + is_block_transposition)

