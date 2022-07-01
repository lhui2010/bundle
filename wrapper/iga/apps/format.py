"""
format related
"""
#!/usr/bin/env python
import re
import argparse
import sys
import collections

from iga.apps.base import emain


class BedIO():
    """ format the following line into bed information of two sides """
#1	1	8468222	8471062	B73_Zm00001d027311_T001	0	+	A188G02321-t1	A188G02321-t1
    # 1-to-1	OrthoFinder+Syntenic	.	=	1	2746424	2749031	B73_Zm00001d027311_T001	.	+
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
                if len(mylist) < 6:
#Some times bed_field is not complete, use '.' to fill them
                    short = 6 - len(mylist)
                    mylist += ["."]*short
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
#            print(k)
#            print(new_key)
            new_bed_line[new_key] = "\t".join([self.chr_start_end[k], new_key, self.score[k], self.strand[k]])
            new_strand[new_key] = self.strand[k]
            new_score[new_key] = self.score[k]
            new_chr_start_end[new_key] = self.chr_start_end[k]
        self.bed_line = new_bed_line
        self.strand = new_strand
        self.score = new_score
        self.chr_start_end = new_chr_start_end

    def print_bed(self):
        """Print entire bed in to a str"""
        print_buf = ""
        for k in self.bed_line:
            print_buf += self.bed_line[k] + "\n"
        return print_buf


def read_ortho(ortho_file):
    this_dict = collections.OrderedDict()
    with open(ortho_file) as fh:
        for line in fh:
            mylist = line.rstrip().split()
            this_dict[mylist[0]] = mylist[1]
    return this_dict


def add_loci_to_ortho(ortho=None, qrybed=None, refbed=None):
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

    ortho_applied = 0

    qry_bed = BedIO(qrybed)
    ref_bed = BedIO(refbed)
    key_list = []

    if(ortho):
        get_ref_by_qry = read_ortho(ortho)
        ortho_applied = 1
        key_list = get_ref_by_qry.keys()
    else:
        key_list = qry_bed.chr_start_end.keys()

#    ref_tobesub = re.compile(r'-R.*')
#    ref_subto = ""
#    ref_bed.rename(ref_tobesub, ref_subto)

    for k in key_list:
        qry_name = k
        assert qry_bed.chr_start_end[qry_name], "No bed_infor for qry {}".format(qry_name)
        loci_qry = qry_bed.chr_start_end[qry_name]
        score = qry_bed.score[qry_name]
        strand_qry = qry_bed.strand[qry_name]

        if(ortho_applied):
            assert get_ref_by_qry[k], "No orthologs for {}".format(k)
            ref_name = get_ref_by_qry[k]
        else:
            ref_name = k
        assert ref_bed.chr_start_end[ref_name], "No bed_infor for ref {}".format(ref_name)
        loci_ref = ref_bed.chr_start_end[ref_name]
        strand_ref = ref_bed.strand[ref_name]

        # Pv01	Phvul.001G003700	211727	214060	+	Aradu.A06	Aradu.DSS3T	17802656	17804120	-
        line = "\t".join([loci_qry, loci_ref, qry_name, score, strand_qry, strand_ref, ref_name])
        print(line)


#Prepare fasta
if __name__ == '__main__':
    emain()
