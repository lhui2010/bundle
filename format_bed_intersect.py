#!/usr/bin/env python

import sys

def intersect_loci(a,b,c,d):
    """
    Loci of intersection 
           start end
    loci1: a     b
    loci2: c     d
    """
    return (max(int(a),int(c)), min(int(b),int(d)))

with open(sys.argv[1]) as fh:
    for line in fh:
        (chra, starta, enda, chrb, startb, endb, overlap_size) = line.rstrip().split()
        tmp = set(range(int(starta), int(enda) + 1)) & set(range(int(startb), int(endb) + 1))
        (new_start, new_end) = intersect_loci(starta, enda, startb, endb)
        #exit()
        print("\t".join([chra, str(new_start), str(new_end)]))
