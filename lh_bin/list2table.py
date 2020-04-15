#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description='Merge different list into a table by keywords')
parser.add_argument('files', metavar='files', type=str, nargs='+',
                    help='filenames need to be merged')

args = parser.parse_args()


total_dict = {}

for this_file in args.files:
    with open (this_file) as fh:
        for line in fh:
            this_list = line.split()
            if(not total_dict.__contains__(this_list[0])):
                total_dict[this_list[0]] = {}
            total_dict[this_list[0]][this_file] = this_list[1]

#Format Print
header = "Key\t" + "\t".join(args.files)
print (header)

for k in (sorted(total_dict.keys())):
    content = k 
    for filename in args.files:
        if(not total_dict[k].__contains__(filename)):
            total_dict[k][filename] = '0'
        content += "\t" + total_dict[k][filename] 
    print(content)

