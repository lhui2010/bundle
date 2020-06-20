#!/usr/bin/env python

import argparse
import textwrap
import logging

import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from assembly_agp2bedpe import AssemblyIO

#Prepare fasta
def main():
    prog_name = "assembly_fasta_to_build"
    usage = "Another python program"

#use of textwrap allowed custome format like "\t\tusage"
    parser = argparse.ArgumentParser(
        prog=prog_name, 
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(usage), 
        epilog="")
    parser.add_argument("assembly", help="final.assembly file")
    parser.add_argument("fasta", help="input fasta for juicer & 3dDNA")
#    parser.add_argument("-f", "--flanking", default=10000, type=int, help="flanking distance default (1000)")
    args = parser.parse_args()  
    assembly_file = args.assembly
    fasta_file = args.fasta
#Read assembly file
    this_assembly_obj = AssemblyIO(assembly_file)
#Read fasta file
    logging.warning("Reading fasta sequences..\n")
    fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    for each_chr in this_assembly_obj.bed_dict:
        chr_seq = ""
        for each_feat in this_assembly_obj.bed_dict[each_chr]:
            chr_seq += str(each_feat.get_fasta(fasta_dict))

        record = SeqRecord(Seq(chr_seq),
                id=each_chr, name='', description='')
        print(record.format('fasta'))
#    print(dir(fasta_dict["Hic.fastq.gz.counts_GATC.20g10"]))
#    print(help(fasta_dict["Hic.fastq.gz.counts_GATC.20g10"].reverse_complement))

if __name__ == "__main__":
    main()
