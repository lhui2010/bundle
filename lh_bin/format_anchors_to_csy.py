#!/usr/bin/env python
import re
import argparse
import sys
import collections



class BedIO():
    """ format the following line into bed information of two sides """
#1    1    8468222    8471062    B73_Zm00001d027311_T001    0    +    A188G02321-t1    A188G02321-t1    1-to-1    OrthoFinder+Syntenic    .    =    1    2746424    2749031    B73_Zm00001d027311_T001    .    +
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
                self.chr_start_end[self.gene] = "\t".join([self.chr, self.gene, self.start, self.end])

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

def read_anchors(anchors_file):
    this_dict = collections.OrderedDict()
    with open(anchors_file) as fh:
        for line in fh:
            mylist = line.rstrip().split()
            this_dict[mylist[0]] = mylist[1]
    return this_dict


def main():
    """ Input example
    First bed:
    ==> NS3.0.A188.bed <==
    1    124014111    124020550    A188G38620-t1    0    +
    1    219433070    219439641    A188G06010-t1    0    -
    1    219690949    219694044    A188G06018-t1    0    +
    1    219847690    219849266    A188G06029-t1    0    -
    ==> NS3.0.B73.bed <==
    5    218008991    218012951    A188G19803-t1-R1-1-A1    149.13    -
    5    217936200    217940096    A188G19806-t1-R1-1-A1    232.826    +
    2    200679459    200679653    A188G19821-t1-R1-1-A1    14.5    -
    4    47177984    47178744    A188G19853-t1-R2-2-A1    23.9494    -
    """

    parser = argparse.ArgumentParser(description='Merge two bed into one bedpe')
    parser.add_argument("QRY_BED", help="Fasta of query fasta files")
    parser.add_argument("REF_BED", help="Fasta of reference fasta files")
    parser.add_argument("-o", "--anchors", default='', help="specifying .anchors file")
    parser.add_argument("-s", "--simple", default='', help="specifying .simple file")
    args = parser.parse_args()  

    QRY_BED = args.QRY_BED
    REF_BED = args.REF_BED
    ANCHORS = args.anchors
    SIMPLE = args.simple
    anchors_applied = 0

    qry_bed = BedIO(QRY_BED)
    ref_bed = BedIO(REF_BED)
    key_list = []

#Pv01 DAGchainer 123 456
    block_prefix = {}
# 54
    block_N_pairs = {}
#matches=chr:start:end
    block_suffix = {}

    program="MCScanX"


    if(SIMPLE):
        simple_dict = {};
        with open(SIMPLE) as fh:
#      Block   Chr     Start   End     Span    StartGene       EndGene GeneSpan        Orientation
#      Cercis_chinensis-Medicago_truncatula-block-0    chr03_Cechi     12363365        12384700        21336   CECHI00021650-t1_Cechi  CECHI00021658-t1
#      Cercis_chinensis-Medicago_truncatula-block-0    MtrunA17CP_Metru        29151   50087   20937   mRNA_MtrunA17CPg0492421_Metru   mRNA_MtrunA17CPg
            fh.readline()
            for line in fh:
                (Block, Chr, Start, End, Span, StartGene, EndGene, GeneSpan, Orientation) = line.rstrip().split()
                if Block in block_prefix:
                    block_suffix[Block] = "matches={}:{}:{}".format(Chr, int(Start)-1, int(End))
                else:
                    block_prefix[Block] = "{}\t{}\t{}\t{}".format(Chr, program, int(Start)-1, int(End))


    if(ANCHORS):
        #### Alignment 0: score=250.0 e_value=5e-11 N=6 MtrunA17Chr0c01_Metru&chr11_Zeins plus
        with open(args.anchors) as fh:
            block_id = -1
            block_prefix_tag = "{}-{}-block-".format(QRY_BED.replace('.bed', ''), REF_BED.replace('.bed', ''))
            for line in fh:
                Orientation = ""
                if(line.startswith("#")):
                    block_id += 1
                    this_block_name = block_prefix_tag + str(block_id)
                    mylist = line.rstrip().split()
                    if(mylist[-1] == 'plus'):
                        Orientation = '+'
                    elif(mylist[-1] == 'minus'):
                        Orientation = '-'
                    else:
                        print("Unknown orientation for {}".format(line))
                        exit(1)
                    gene_pair_count = mylist[-3].replace('N=', '')
                    print("{}\t{}\t{}\t{}".format(block_prefix[this_block_name], gene_pair_count, Orientation, block_suffix[this_block_name]) + ";median_Ks=6")
                else:
                    (qry_gene, ref_gene, score) = line.rstrip().split()
                    new_line = "{}\t{}".format(qry_bed.chr_start_end[qry_gene], ref_bed.chr_start_end[ref_gene])
                    print(new_line)

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

