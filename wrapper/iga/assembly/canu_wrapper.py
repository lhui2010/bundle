#!/usr/bin/env python

import argparse
import textwrap

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
    with open(qry1_file) as fh:
        for line in fh:
            mylist = line.rstrip().split()

if __name__ == "__main__":
    main()
