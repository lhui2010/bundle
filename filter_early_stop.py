#Zm00001d013095_T001|000154F|arrow|arrow|pilon:1448735..1452656|-|gene cover:1065(100%)|score:532.754|rank:1
#HSP_ID[6]:(1452571-1452656);query:(1-86); pid: 100
#HSP_ID[4]:(1452369-1452481);query:(85-197); pid: 100
#HSP_ID[11]:(1452217-1452251);query:(194-228); pid: 100
#HSP_ID[2]:(1451704-1451872);query:(229-397); pid: 100
#HSP_ID[1]:(1451105-1451288);query:(396-579); pid: 100
#
#Zm00001d009939_T001|001704F|arrow|arrow|pilon:30526..31453|+|gene cover:618(100%)|score:308.975|rank
#HSP_ID[2]:(31217-31453);query:(382-618); pid: 100
#HSP_ID[3]:(30895-31050);query:(237-392); pid: 99.359
#HSP_ID[1]:(30526-30763);query:(1-238); pid: 100

import sys
import re
from Bio import SeqIO

#xx.py genblast genome.fasta

#TODO replace sys with parse arg

import argparse

def main():
    parser = argparse.ArgumentParser(description='Extract cds, bed, and gff from genblasta result')
    parser.add_argument('FA', type=str, nargs = 1,
                        help='Genome.fasta')
    parser.add_argument('--version', action='version', version='%(prog)s 0.1')

    args = parser.parse_args()

    ffasta=args.FA[0]

    raw_dict = SeqIO.to_dict(SeqIO.parse(ffasta, "fasta"))            
    genome_dict = {}
#rename genome_dict
    for old_key in raw_dict.keys():
        if(re.search(r'\*|U', raw_dict[old_key].seq.__str__())):
            continue
        else:
            genome_dict[old_key] = raw_dict[old_key]

    with open(ffasta + "filterStop.fa", "w") as output_handle:
        for k in genome_dict.keys():
            SeqIO.write(genome_dict[k], output_handle, "fasta")

main()
