#!/usr/bin/env python
import re

#python ${SCRIPT_DIR}/lift_gff.py ../${DIFF_TABLE}.${i} syn_genblast.gff.raw > syn_genblast.gff

class GFFIO():
    """simple gff class"""
    def __init__(self, filename):
        with open(filename) as fh:
            self.GFF_list = fh.readlines()

    def liftover_by_chrid(self):
        for line_id,line in enumerate(self.GFF_list):
#1_265249172_266462661   genBlastG       transcript      1203511 1206816 55.1412 +       .       ID=A188_A188G
            mylist = line.rstrip().split()
            (newchrid, cut_start, cut_end) = mylist[0].split('_')
            mylist[0] = newchrid
            mylist[3] = str(int(mylist[3]) + int(cut_start) - 1)
            mylist[4] = str(int(mylist[4]) + int(cut_start) - 1)
            self.GFF_list[line_id] = "\t".join(mylist) + "\n"
    def print_gff(self):
        for line in self.GFF_list:
            print(line.rstrip())

#def read_offset_file(filename, offset, chrid):
#    with open(filename) as fh:
#        for line in fh:
#            (name, chrid, start, end) = line.rstrip().split()
#            offset[name] = start
#            chrid[name] = chrid

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Lift-over gff coordinates from genblast gff to genome loci.')
#    parser.add_argument("offset", metavar="<chr1_gene_diff.tab.out.A188>", type=str, help="file describing where genblast queries were cut")
    parser.add_argument("gff", metavar="<genes.gff>", type=str, help="Gff file to be lifted-over")
    args = parser.parse_args()
    gff_obj = GFFIO(args.gff)
    gff_obj.liftover_by_chrid()
    gff_obj.print_gff()


