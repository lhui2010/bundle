#!/usr/bin/env python
import re
import argparse
import sys

#Args
parser = argparse.ArgumentParser(description='Sort 2nd file by order of first file')
parser.add_argument('file_key', type=str, nargs = 1,
                    help='file of key')
parser.add_argument('file_table', type=str, nargs = 1,
                    help='file of table')
parser.add_argument('--version', action='version', version='%(prog)s 0.1')
args = parser.parse_args()

#read 2nd into a dict
table_dict = dict()
with open (args.file_table[0]) as fh:
    for line in fh:
        mylist = line.split()
        table_dict[mylist[0]] = line.rstrip()

order=[]
with open (args.file_key[0]) as fh:
    for line in fh:
        mykey = line.rstrip()
        if(table_dict.__contains__(mykey)):
            print(table_dict[mykey])
        order.append(line.rstrip())
        
