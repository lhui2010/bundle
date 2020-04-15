#!/usr/bin/env python
import argparse
#import sys

#Args
parser = argparse.ArgumentParser(description='python template')
parser.add_argument('INPUT', type=str, nargs = 1,
                    help='input file')
parser.add_argument('--version', action='version', version='%(prog)s 0.1')
args = parser.parse_args()

inputs = args.INPUT[0]
