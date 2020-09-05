#!/usr/bin/env python

import argparse
import textwrap
import re
import os

def main():
    prog_name = "Another python program"
    usage = """
    Input:
    ==> 140309at3193.gff <==
    000017F|arrow_np1212    AUGUSTUS        gene    4073034 4074064 0.28    -       .       g1
    000017F|arrow_np1212    AUGUSTUS        transcript      4073034 4074064 0.28    -       .       g1.t1
    000017F|arrow_np1212    AUGUSTUS        intron  4073034 4073049 0.82    -       .       transcript_id   "g1.t1";        gene_id "g1";
    000017F|arrow_np1212    AUGUSTUS        intron  4073078 4073169 0.98    -       .       transcript_id   "g1.t1";        gene_id "g1";
    000017F|arrow_np1212    AUGUSTUS        intron  4073290 4073922 0.81    -       .       transcript_id   "g1.t1";        gene_id "g1";
    000017F|arrow_np1212    AUGUSTUS        CDS     4073050 4073077 0.82    -       2       transcript_id   "g1.t1";        gene_id "g1";
    000017F|arrow_np1212    AUGUSTUS        exon    4073050 4073077 .       -       .       transcript_id   "g1.t1";        gene_id "g1";
    000017F|arrow_np1212    AUGUSTUS        CDS     4073170 4073289 0.82    -       2       transcript_id   "g1.t1";        gene_id "g1";
    000017F|arrow_np1212    AUGUSTUS        exon    4073170 4073289 .       -       .       transcript_id   "g1.t1";        gene_id "g1";
    000017F|arrow_np1212    AUGUSTUS        CDS     4073923 4074022 0.71    -       0       transcript_id   "g1.t1";        gene_id "g1";

    ==> 159907at3193.gff <==
    000027F|arrow_np1212    AUGUSTUS        gene    1027988 1032880 0.01    -       .       g1
    000027F|arrow_np1212    AUGUSTUS        transcript      1027988 1032880 0.01    -       .       g1.t1
    000027F|arrow_np1212    AUGUSTUS        tts     1027988 1027988 .       -       .       transcript_id   "g1.t1";        gene_id "g1";
    """

    parser = argparse.ArgumentParser(
        prog=prog_name, 
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(usage), 
        epilog="")
    parser.add_argument("qry1", help="qry1 file", nargs='+')
#    parser.add_argument("-f", "--flanking", default=10000, type=int, help="flanking distance default (1000)")
    args = parser.parse_args()  

    qry1_file = args.qry1
#    flanking_distance = args.flanking
    for q in qry1_file:
        with open(q) as fh:
            bsname=os.path.splitext(q)[0]
            count = 1
            tt_name = ''
            for line in fh:
                mylist = line.rstrip().split(maxsplit=8)
#    000017F|arrow_np1212    AUGUSTUS        transcript      4073034 4074064 0.28    -       .       g1.t1
#    000017F|arrow_np1212    AUGUSTUS        intron  4073034 4073049 0.82    -       .       transcript_id   "g1.t1";        gene_id "g1";
                if(mylist[2] == "transcript"):
                    mylist[2] = re.sub('transcript', 'match', mylist[2])
                    count = 1
                    tt_name = bsname + mylist[-1]
                    mylist[-1] = "ID=" + tt_name + ";" + "Name=" + tt_name
                elif(mylist[2] == "CDS"):
                    mylist[2] = re.sub('CDS', 'match_part', mylist[2])
                    mylist[-1] = "ID=" + tt_name + ":" +str(count) + ";Parent=" + tt_name
                    count += 1
                else:
                    continue
                print("\t".join(mylist))

if __name__ == "__main__":
    main()
