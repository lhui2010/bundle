import sys
import os
import argparse
#from get_unsyntenic_genes_syntenic_loci import Bedpe

class Bedpe_original():
    """ 
    format the following line into bed information of two sides 
    1	124014111	124020550	1	120051289	120060421	A188G38620-t1	0	+	-
    1	219433070	219439641	9	152298055	152298468	A188G06010-t1	0	-	-
    """
    def __init__(self, line):
        """Initialize the values"""
        mylist = line.rstrip('\n').split('\t')
        #print(mylist)
        self.gene_left = mylist[6]
        (self.chr_left, self.start_left, self.end_left) = mylist[0:3]
        self.strand_left = mylist[8]
        if(len(mylist) > 10):
            self.gene_right = mylist[10]
        else:
            self.gene_right = mylist[6]
        (self.chr_right, self.start_right, self.end_right) = mylist[3:6]
        self.strand_right = mylist[9]
        self.score = mylist[7]

    def get_left_bed(self, flanking=0):
        new_start = str(int(self.start_left) - flanking)
        new_end = str(int(self.end_left) + flanking)
        line = "\t".join([self.chr_left, new_start, 
        new_end, self.gene_left, self.score, self.strand_left])
        return line

    def get_right_bed(self, flanking=0):
        new_start = str(int(self.start_right) - flanking)
        new_end = str(int(self.end_right) + flanking)
        line = "\t".join([self.chr_right, new_start, 
        new_end, self.gene_right, self.score, self.strand_right])
        return line

    def get_left_loci(self, flanking=0):
        new_start = str(int(self.start_left) - flanking)
        new_end = str(int(self.end_left) + flanking)
        line = self.chr_left + ":" + new_start + "-" + new_end
        return line

    def get_right_loci(self, flanking=0):
        new_start = str(int(self.start_right) - flanking)
        new_end = str(int(self.end_right) + flanking)
        line = self.chr_right + ":" + new_start + "-" + new_end 
        return line


def main():
    usage = """Plot genome track from bedpe file"""

    parser = argparse.ArgumentParser(usage)
    parser.add_argument("bed_pe", help="bed pe file")
    parser.add_argument("bed_left", help="bed file for left gene list")
    parser.add_argument("bed_right", help="bed file for right gene list")
    parser.add_argument("-f", "--flanking", default=10000, type=int, help="flanking distance default (1000)")
    args = parser.parse_args()  

    bedpe_file = args.bed_pe
    bedleft_file = args.bed_left
    bedright_file = args.bed_right
    flanking_distance = args.flanking

    file_bed_total = {}
    file_bed_total["left"] = bedleft_file
    file_bed_total["right"] = bedright_file
    with open(bedpe_file) as fh:
        for line in fh:
            this_bedpe = Bedpe_original(line)
            this_gene = this_bedpe.gene_left
            this_dict = {}
            this_dict["leftloci"] = this_bedpe.get_left_loci(flanking_distance)
            this_dict["leftbed"] = this_bedpe.get_left_bed(flanking_distance)
            this_dict["rightloci"] = this_bedpe.get_right_loci(flanking_distance)
            this_dict["rightbed"] = this_bedpe.get_right_bed(flanking_distance)
#            print(this_bedpe.get_left_loci())
#            print(this_bedpe.get_right_loci())

            for side in ["left", "right"]:
                output = this_gene + "." + side
                bed_content = this_dict[side + "bed"]
                loci_content = this_dict[side + "loci"]
                trackname = output + ".ini"
                file_loci_bed = output + "loci.bed"
                file_bed = output + ".bed"
                file_ini = output + ".ini"
                with open(output + "loci.bed", 'w') as fh:
                    fh.write(bed_content)
                with open(output + ".ini", 'w') as fh:
                    ini_file = """
[genes]
file = """ + file_bed + """
height = 7
title = genes
fontsize = 10
file_type = bed
gene_rows = 10
color = Reds

[x-axis]
fontsize=10
"""
                    fh.write(ini_file)                
                os.system("bedtools intersect -a {} -b {} -wb |cut -f7,8,9,10,11,12 |bedtools sort -i - > {}".format(file_loci_bed, file_bed_total[side], file_bed))
#                os.system("make_tracks_file --trackFiles {} -o {}".format("test.bed" ,trackname))
                os.system("pyGenomeTracks --trackLabelFraction 0.001 --tracks {} -o {} --region {}".format(trackname, output + ".pdf", loci_content))
#            break

if __name__ == '__main__':
    main()
