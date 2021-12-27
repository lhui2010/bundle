#!/usr/bin/env python
import re
import argparse
import sys


class BedIO():
    """ format the following line into bed information of two sides """
#1	1	8468222	8471062	B73_Zm00001d027311_T001	0	+	A188G02321-t1	A188G02321-t1	1-to-1	OrthoFinder+Syntenic	.	=	1	2746424	2749031	B73_Zm00001d027311_T001	.	+
    def __init__(self, filebed):
        """Initialize the values"""
        self.bed_line = {}
        self.strand = {}
        self.score = {}
        self.chr_start_end = {}
        """Read filebed into a bed dict"""
        with open(filebed) as fh:
            for line in fh:
                mylist = line.rstrip('\n').split('\t')
        #print(mylist)
                (self.chr, self.start, self.end, self.gene, self.this_score, self.this_strand) = mylist[0:6]
                self.bed_line[self.gene] = "\t".join(mylist[0:6])
                self.strand[self.gene] = self.this_strand
                self.score[self.gene] = self.this_score
                self.chr_start_end[self.gene] = "\t".join(mylist[0:3])

    def rename(self, re_tobesub, re_subto):
        """Change the keys with re.sub(re_tobesub, re_subto)"""
        new_bed_line = {}
        new_strand = {}
        new_score = {}
        new_chr_start_end = {}
        for k in self.bed_line:
            new_key = re.sub(re_tobesub, re_subto, k)
            new_bed_line[new_key] = self.bed_line[k]
            new_strand[new_key] = self.strand[k]
            new_score[new_key] = self.score[k]
            new_chr_start_end[new_key] = self.chr_start_end[k]
        self.bed_line = new_bed_line
        self.strand = new_strand
        self.score = new_score
        self.chr_start_end = new_chr_start_end

def main():
    """ Input example
    First bed:
    ==> NS3.0.A188.bed <==
    1	124014111	124020550	A188G38620-t1	0	+
    1	219433070	219439641	A188G06010-t1	0	-
    1	219690949	219694044	A188G06018-t1	0	+
    1	219847690	219849266	A188G06029-t1	0	-
    ==> NS3.0.B73.bed <==
    5	218008991	218012951	A188G19803-t1-R1-1-A1	149.13	-
    5	217936200	217940096	A188G19806-t1-R1-1-A1	232.826	+
    2	200679459	200679653	A188G19821-t1-R1-1-A1	14.5	-
    4	47177984	47178744	A188G19853-t1-R2-2-A1	23.9494	-
    """

    parser = argparse.ArgumentParser(description='Regex based exact sequence alignment tool\nResult is 1-based')
    parser.add_argument('QRYBED', type=str, nargs = 1,
                        help='BED of query')
    parser.add_argument('REFBED', type=str, nargs = 1,
                        help='BED of reference')
    parser.add_argument('--version', action='version', version='%(prog)s 0.1')
    args = parser.parse_args()

    qry_bed = BedIO(args.QRYBED[0])
    ref_bed = BedIO(args.REFBED[0])

    ref_tobesub = re.compile(r'-R.*')
    ref_subto = ""

    ref_bed.rename(ref_tobesub, ref_subto)

    for k in qry_bed.chr_start_end.keys():
        loci_qry = qry_bed.chr_start_end[k]
        loci_ref = ref_bed.chr_start_end[k]
        name = k
        score = qry_bed.score[k]
        strand_qry = qry_bed.strand[k]
        strand_ref = ref_bed.strand[k]
        line = "\t".join([loci_qry, loci_ref, name, score, strand_qry, strand_ref])
        print(line)


#Prepare fasta
if __name__ == '__main__':
    main()

#@sys.stderr.write("Reading fasta sequences..\n")
#@with open(args.QRYBED[0]) as fh:
#@    for line in fh:
#@        mylist = line.rstrip().split()
#@qry_dict = SeqIO.to_dict(SeqIO.parse(args.QRYBED[0], "fasta"))
#@ref_dict = SeqIO.to_dict(SeqIO.parse(args.REFBED[0], "fasta"))
#@ref_dict_reverse = {}
#@print_buf = ""
#@
#@#Loop
#@ref_fa_len = {}
#@sys.stderr.write("Reversing Ref fasta sequences..\n")
