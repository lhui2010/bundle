#!/usr/bin/env python
import re
#import regex as re
import argparse
from Bio import SeqIO
import sys
#import pprint

#Args
parser = argparse.ArgumentParser(description='Find Gap (nN) in a fasta file and output their position in bed format (start is 0-based)')
parser.add_argument('QRY', type=str, nargs = 1,
                    help='fasta of query sequences')
parser.add_argument('--version', action='version', version='%(prog)s 0.1')
args = parser.parse_args()

#Prepare fasta
sys.stderr.write("Reading fasta sequences..\n")
qry_dict = SeqIO.to_dict(SeqIO.parse(args.QRY[0], "fasta"))
print_buf = ""

#Loop
sys.stderr.write("Reversing Ref fasta sequences..\n")

key = re.compile(r'[nN]+')

for qry_id in qry_dict.keys():
    result_line = ""
    sys.stderr.write("Searching " + str(key) + " in " + qry_id + "..\n")
    for result in key.finditer(qry_dict[qry_id].seq.__str__()):
#    result=re.findall(r'[nN]+', qry_dict[qry_id].seq.__str__(), re.IGNORECASE)
#        print(result.start(), result.group())
        ncount = len(result.group())
        gap_ID = "N{" + str(ncount) + "}"
        if(result != None):
            strand = '+'
            result_line= "\t".join([qry_id,str(result.start()),
                str(result.end()), gap_ID ,".",strand])
        if(result_line != ""):
            print_buf += result_line + "\n"
print(print_buf, end ='')

