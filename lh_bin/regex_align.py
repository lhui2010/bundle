#!/usr/bin/env python3
import re
import argparse
from Bio import SeqIO
import sys

#Args
parser = argparse.ArgumentParser(description='Regex based exact sequence alignment tool\nResult is 1-based')
parser.add_argument('QRY', type=str, nargs = 1,
                    help='fasta of query sequences')
parser.add_argument('REF', type=str, nargs = 1,
                    help='fasta of reference sequences')
parser.add_argument('--version', action='version', version='%(prog)s 0.1')
args = parser.parse_args()

#Prepare fasta
sys.stderr.write("Reading fasta sequences..\n")
qry_dict = SeqIO.to_dict(SeqIO.parse(args.QRY[0], "fasta"))
ref_dict = SeqIO.to_dict(SeqIO.parse(args.REF[0], "fasta"))
ref_dict_reverse = {}
print_buf = ""

#Loop
ref_fa_len = {}
sys.stderr.write("Reversing Ref fasta sequences..\n")
for ref_id in ref_dict:
    ref_dict_reverse[ref_id] = ref_dict[ref_id].seq.reverse_complement()
    ref_fa_len[ref_id] = len(ref_dict[ref_id].seq)

for qry_id in qry_dict.keys():
    result_line = ""
    for ref_id in ref_dict.keys():
        sys.stderr.write("Searching " + qry_id + " in " + ref_id + "..\n")
        result=re.search(qry_dict[qry_id].seq.__str__(), ref_dict[ref_id].seq.__str__(), re.IGNORECASE)
        if(result != None):
            strand = '+'
            result_line= "\t".join([ref_id,str(result.start()+1),
                str(result.end()), qry_id,".",strand])
        else:
            result=re.search(qry_dict[qry_id].seq.__str__(), 
                ref_dict_reverse[ref_id].__str__(),re.IGNORECASE)
            if(result != None):
                strand = '-'
                rev_start = ref_fa_len[ref_id] - result.start()
                rev_end = ref_fa_len[ref_id] - result.end() + 1
                result_line= "\t".join([ref_id,str(rev_end),
                    str(rev_start), qry_id,".",strand])
        if(result_line != ""):
            print_buf += result_line + "\n"
print(print_buf, end ='')

