#!/usr/bin/env python

import argparse
import textwrap
import re


def main():
    prog_name = "Another python program"
    usage = "Another python program"

    parser = argparse.ArgumentParser(
        prog=prog_name,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(usage),
        epilog="")
    parser.add_argument("qry1", help="qry1 file")
    parser.add_argument("qry2", help="qry2 file")
    #    parser.add_argument("-f", "--flanking", default=10000, type=int, help="flanking distance default (1000)")
    args = parser.parse_args()

    qry1_file = args.qry1
    qry2_file = args.qry2
    #    flanking_distance = args.flanking

    keyword_dict = {}
    with open(qry1_file) as fh:
        for line in fh:
            mylist = line.rstrip().split()
            keyword_dict[mylist[0]] = 1

    with open(qry2_file) as fh:
        for line in fh:
            flag = 0
            my_gene = ''
            mylist = line.rstrip().split('\t', 8)
            if('Parent' in mylist[-1]):
                my_gene = re.sub(r'.*Parent=', '', mylist[-1])
                my_gene = re.sub(r';.*', '', my_gene)
            elif('ID' in mylist[-1]):
                my_gene = re.sub(r'.*ID=', '', mylist[-1])
                my_gene = re.sub(r';.*', '', my_gene)
            if(my_gene in keyword_dict):
                print(line.rstrip())

if __name__ == "__main__":
    main()
