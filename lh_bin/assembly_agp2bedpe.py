#!/usr/bin/env python

#@import sys
#@# insert at 1, 0 is the script path (or '' in REPL)
#@sys.path.insert(1, '/ds3200_1/users_root/yitingshuang/lh/bin/bundle/curate_orthogene/scripts')
#@from add_score import teBedIO

import argparse
import re
import textwrap
from pprint import pprint

class BedFeat():
    """
    An object to store each bed line
    ==var==     ==type==
    chromosome  str
    start       int
    end         int
    name        str
    strand      str
    score       str
    """
    def __init__(self, assembly_line='', chromosome='.', start='0',\
    end = '0', name='.', strand = '.', score = '.'):
        """ 
        init seqfrag object with .0.review.assembly line 
        (1-based in both start and end loci)
        """
        (self.chromosome, self.start, self.end, self.name, self.strand, self.score) \
            = [chromosome, int(start), int(end), name, strand, score]
#        if assembly_line != '':
#           myl = assembly_line.rstrip().split()
#           (name, start, end) = myl
#        self.chromosome = name.replace('>', '')
        self.start = int(start)
        self.end = int(end)

    def get_line(self):
        return "\t".join([self.chromosome, str(self.start), str(self.end), self.name, self.score, self.strand])

    def get_fasta(self, fa_dict):
        """
        return a fasta of a segment like bedtools getfasta, 
        where self is a bed line and fa_dict is a dict of 
        Bio.SeqRecord
        """
#zero based coordinate
        zero_start = self.start - 1
        zero_end = self.end -1 + 1
        fa_seq = fa_dict[self.chromosome][zero_start:zero_end]
        if self.strand == '-':
            fa_seq = fa_seq.reverse_complement()
        return str(fa_seq.seq)

class AGPIO():
    def __init__(self, fileagp):
        """
        Read AGP file into BedIO object
        ==> groups.agp <==
        Hic.fastq.gz.counts_GATC.20g1	1	79244	1	W	tig00005006|arrow_np1212	1	79244	+
        Hic.fastq.gz.counts_GATC.20g1	79245	79344	2	U	100	contig	yes	map
        Hic.fastq.gz.counts_GATC.20g1	79345	123128	3	W	tig00005007|arrow_np1212	1	43784	-
        """
        self.bed_dict = {}
        with open(fileagp) as fh:
            for line in fh:
                (this_chr, this_start, this_end, this_name, this_strand) = [""]*5
                mylist = line.rstrip('\n').split()
                this_score = '.'
                if mylist[4] == "U":
                    (this_chr, this_start, this_end  ) = mylist[0:3]
                    this_name = "GAP100"
                    this_strand = '+'
                else:
                    (this_chr, this_start, this_end, order, seq_type,  this_name, \
                        rel_start, rel_end, this_strand) = mylist
                """Ignore lines with empty name"""
                newfeat = BedFeat(chromosome=this_chr, start=this_start, end=this_end, \
                    name=this_name, strand=this_strand)

                if this_chr not in self.bed_dict:
                    self.bed_dict[this_chr] = []
                self.bed_dict[this_chr].append(newfeat)

    def get_feats_from_range(self, bedfeat):
        """
        Given chr1:1-200 and provided with chr1:2-30 contig1 chr1:31-50 contig2
        Return a list of BedFeat objects
        """
        reverse_strand = {"-":"+", "+":"-"}
        new_name_tag = ""
        feat_list = []
        seperator = "::"
        if bedfeat.chromosome in self.bed_dict:
            for loop in self.bed_dict[bedfeat.chromosome]:
#If no truncating of contig need 
                bed_len = loop.end - loop.start + 1
                new_start = 1
                new_end = new_start + bed_len - 1
                #print("0::" + "\t".join([bedfeat.chromosome, str(bedfeat.start), 
#                    str(bedfeat.end), loop.name, str(loop.start), str(loop.end)]))
                if loop.start < bedfeat.start and loop.end >= bedfeat.start:
                    # bedfeat   =====
                    # loop1   ****
                    # loop2   *********
                    new_start = bedfeat.start - loop.start + 1
                    if loop.end > bedfeat.end:
                        new_end = bedfeat.end - loop.start + 1
                    else:
                        new_end = loop.end - loop.start + 1
                    #print("1::" + "\t".join([loop.name, str(loop.start), str(loop.end)]))
                elif loop.start >= bedfeat.start and loop.start <= bedfeat.end:
                    # bedfeat  =======
                    # loop1      ****
                    # loop2      *********
                    if loop.end > bedfeat.end:
                        new_end = bedfeat.end - loop.start + 1
                    #print("2::" + "\t".join([loop.name, str(loop.start), str(loop.end)]))
                else:
                    continue
                if bedfeat.strand == '+':
                    new_strand = loop.strand
                elif bedfeat.strand == '-':
                    new_strand = reverse_strand[loop.strand]
                #ctg01::100::1::80 contig_name::original_contig_len::contig_start::contig_end
                new_name_tag = seperator.join([loop.name, str(bed_len), str(new_start), \
                    str(new_end)])
                new_feat = BedFeat(start='1', end=str(new_end - new_start + 1), strand=new_strand, \
                    name=new_name_tag)
                feat_list.append(new_feat)
            return feat_list
        else:
            return "ERROR: %s not found in AGP file" % bedfeat.chromosome

class AssemblyIO():
    def __init__(self, assembly_file):
        """
        Input assembly file
        #==> allhic.0.review.assembly <==
        #>Hic.fastq.gz.counts_GATC.20g1:::fragment_1 1 5423057
        #>Hic.fastq.gz.counts_GATC.20g1:::fragment_2:::debris 2 50000
        #>Hic.fastq.gz.counts_GATC.20g1:::fragment_3 3 7785000
        #...
        #1 -3 13 22 -5
        """
        self.bed_dict = {}
        initial_chr_id = 1
        with open(assembly_file) as fh:
#Store contig name before hic
            seg_list = [""]
#Store groups assemblied by JCBAT
            chr_list = [""]
            chr_id = initial_chr_id
            last_chr_id = ''
            offset = 0

            for line in fh:
                (this_chr, this_start, this_end, this_name, this_strand) = ["."]*5
#Each line is a new chromosome
                if line.startswith(">"):
                    (frag_name, frag_order, frag_len) = line.rstrip().split()
                    frag_name = frag_name[1:]
                    frag_len = int(frag_len)
                    if 'fragment' in frag_name:
                        allhic_chr_name = re.sub(r':.*', '', frag_name)
                    else:
                        allhic_chr_name = frag_name
                    if allhic_chr_name != last_chr_id:
                        offset = 0
                    new_start = offset + 1
                    new_end =  offset + frag_len
                    seg_list.append(BedFeat(chromosome=allhic_chr_name, start=new_start,
                        end = new_end))
                    last_chr_id = allhic_chr_name
                    offset += frag_len
                else:
                    this_chr = "chr" + "%02d" % chr_id
                    #print(this_chr)
                    chr_id += 1
                    if this_chr not in self.bed_dict:
                        self.bed_dict[this_chr] = []
#1 43 -21
                    mylist = line.rstrip().split()
                    for i in mylist:
                        this_strand = '+'
                        this_id=int(i)
                        if this_id < 0:
                            this_strand = '-'
                            this_id = int(i[1:])
                        translated_seg = seg_list[this_id]
                        translated_seg.strand = this_strand
                        self.bed_dict[this_chr].append(translated_seg)
        
def main():
    usage = """
    Another python program
    Input:
        ==> allhic.0.review.assembly <==
        >Hic.fastq.gz.counts_GATC.20g1:::fragment_1 1 5423057
        >Hic.fastq.gz.counts_GATC.20g1:::fragment_2:::debris 2 50000
        >Hic.fastq.gz.counts_GATC.20g1:::fragment_3 3 7785000
        ...
        1 -3 13 22 -5
        ==> groups.agp <==
        Hic.fastq.gz.counts_GATC.20g1	1	79244	1	W	tig00005006|arrow_np1212	1	79244	+
        Hic.fastq.gz.counts_GATC.20g1	79245	79344	2	U	100	contig	yes	map
        Hic.fastq.gz.counts_GATC.20g1	79345	123128	3	W	tig00005007|arrow_np1212	1	43784	-
    Output:
        NewChr1 1    3000 contig::1_3000 . +
    Example:
        python assembly_agp2bedpe.py allhic.0.review.assembly allhic.groups.agp >CANUcontig_to_JCBATchromsome.bed
    """
    parser = argparse.ArgumentParser(
        prog="assembly_agp_to_bed", 
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(usage), 
        epilog="")
    parser.add_argument("assembly", help="assembly file")
    parser.add_argument("AGP", help="AGP file")
#    parser.add_argument("-f", "--flanking", default=10000, type=int, help="flanking distance default (1000)")
#    flanking_distance = args.flanking
    args = parser.parse_args()  
    assembly_file = args.assembly
    AGP_file = args.AGP

#FileIO
    agp_obj = AGPIO(AGP_file)
    assembly_obj = AssemblyIO(assembly_file)
    new_assembly_dict = {}

#Loop each chromosome
    for chrid_in_assembly in sorted(assembly_obj.bed_dict.keys()):
        offset_new_chr = 0
#Loop each range
        for chr_range in assembly_obj.bed_dict[chrid_in_assembly]:
            new_assembly_dict[chrid_in_assembly] = []
#            pprint(vars(chr_range))
            print("\t".join([chrid_in_assembly, str(chr_range.start), str(chr_range.end), chr_range.chromosome]))
            feats_in_range = agp_obj.get_feats_from_range(chr_range)
#            print(feats_in_range)
            for translated_feat in feats_in_range:
#                print(translated_feat.chromosome)
#                print(chrid_in_assembly)
                translated_feat.chromosome = chrid_in_assembly
                translated_feat.start += offset_new_chr
                translated_feat.end += offset_new_chr
                offset_new_chr += translated_feat.end - translated_feat.start + 1
                new_assembly_dict[chrid_in_assembly].append(translated_feat)
#Debug info                if translated_feat.start == 1:
#@                    print("\t\t" + "\t".join([chr_range.get_line(), translated_feat.chromosome, str(translated_feat.start), str(offset_new_chr)]))
                print(translated_feat.get_line())

if __name__ == "__main__":
    main()
