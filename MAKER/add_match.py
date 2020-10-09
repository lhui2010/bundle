#!/usr/bin/env python

import argparse
import textwrap
import re
from parse import *
import logging

def main():
    prog_name = "Another python program"
    usage = "Another python program"
    """
In pipelines:
#for i in L*.fa; do bsub512 "minimap2 -t20 -C5 -ax splice CORNE.contig.fa ${i} |samtools view -b - | bedtools bamtobed -split -i - > ${i}.bed"; done
#conda activate EDTA
#for i in *bed; do bsub -q Q104C512G_X4  -o output.${i%.bed} -e error.${i%.bed} -J ${i%.bed} "gt bed_to_gff3 ${i} |sort -k9,9 -k1,1 -k7,7 -k4,4n > ${i}.rawgff && python add_match.py ${i} > ${i%.bed}.gff 2> ${i%.bed}.intron_error " ; done
#for i in *bed; do bsub -q Q104C512G_X4  -o output.${i%.bed} -e error.${i%.bed} -J ${i%.bed} "python add_match.py ${i}.rawgff > ${i%.bed}.gff 2> ${i%.bed}.intron_error " ; done
#sort -k9,9 -k1,1 -k7,7 -k4,4n tmp.out > tmp.out.sort
#python add_match.py tmp.out.sort > tmp.out.gff 2> tmp.out.err

Input Example
#000028F|arrow_np1212	.	BED_feature	25928	26136	0	-	.	Name=TRINITY_GG_1_c0_g1_i1
#000053F|arrow_np1212	.	BED_feature	833531	833739	0	-	.	Name=TRINITY_GG_1_c0_g1_i1
#000093F|arrow_np1212	.	BED_feature	21443	21651	0	-	.	Name=TRINITY_GG_1_c0_g1_i1
#000093F|arrow_np1212	.	BED_feature	112492	112699	0	+	.	Name=TRINITY_GG_1_c0_g1_i1
#000188F|arrow_np1212	.	BED_feature	8517	8729	0	-	.	Name=TRINITY_GG_2_c0_g1_i3
#000062F|arrow_np1212	.	BED_feature	27655	27867	0	+	.	Name=TRINITY_GG_2_c0_g1_i3
#000035F|arrow_np1212	.	BED_feature	55447	55659	0	-	.	Name=TRINITY_GG_2_c0_g1_i3
#000035F|arrow_np1212	.	BED_feature	28344	28556	0	+	.	Name=TRINITY_GG_2_c0_g1_i3
#000108F|arrow_np1212	.	BED_feature	28394	28606	0	+	.	Name=TRINITY_GG_2_c0_g1_i3
#000108F|arrow_np1212	.	BED_feature	85313	85524	0	-	.	Name=TRINITY_GG_2_c0_g1_i3
    """

    parser = argparse.ArgumentParser(
        prog=prog_name, 
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(usage), 
        epilog="")
    parser.add_argument("qry1", help="qry1 file")
#    parser.add_argument("-f", "--flanking", default=10000, type=int, help="flanking distance default (1000)")
    args = parser.parse_args()  

    qry1_file = args.qry1
    prefix = re.sub(r'\..*', '', qry1_file)
    intron_cutoff = 20000
    source = 'Trinity_Minimap'
#    flanking_distance = args.flanking
    last_chr = ''
    last_end = ''
    last_strand = ''
    last_pos = []
    last_lines = []
    last_name = ''
    count = 0
    with open(qry1_file) as fh:
        for line in fh:
            if(line.startswith('#')):
                continue
            mylist = line.rstrip().split()
            #print(mylist[-1])
#Name feats
            this_chr = mylist[0]
            this_type = 'match_part'
            this_start = mylist[3]
            this_end = mylist[4]
            this_score = mylist[5]
            this_strand = mylist[6]
            this_phase = mylist[7]
            feat_parse = parse("Name={name}", mylist[-1])
            try:
                feat_name = prefix + feat_parse['name'] + this_chr.replace('|', '') + this_start
            except TypeError:
                logging.warning("TypeError on feat {}".format(mylist[-1]))
            this_feat = "ID={}-exon-{};Name={};Parent={}".format(feat_name, str(count), feat_name, feat_name)
            line = "\t".join([this_chr, this_type, source, this_start, this_end, this_score, this_strand, this_phase, this_feat])

            if(last_chr == '' or 
                    last_name == feat_name and 
                    last_chr == this_chr and 
                    last_strand == this_strand and 
                    abs(last_end - int(this_start) < intron_cutoff)):
                count += 1

                this_feat = "ID={}-exon-{};Name={};Parent={}".format(feat_name, str(count), feat_name, feat_name)
                line = "\t".join([this_chr, this_type, source, this_start, this_end, this_score, this_strand, this_phase, this_feat])
                last_chr = this_chr
                last_end = int(this_end)
                last_pos.append(int(this_start))
                last_pos.append(int(this_end))
                last_strand = this_strand
                last_lines.append(line)
                last_name = feat_name

            else:
                if(last_name == feat_name and last_chr == this_chr and last_strand == this_strand):
                    logging.warning("LargeIntron {} on {}".format(str(abs(last_end - int(this_start))), feat_name ))
#Prepare print
                match_chr = last_chr
                match_type = "match"
                match_start = str(min(last_pos))
                match_end = str(max(last_pos))
                match_score = '0'
                match_strand = last_strand
                match_phase = '.'
                match_feat = "ID={};Name={}".format(last_name, last_name)
                match_line = "\t".join([match_chr, match_type, source, match_start, match_end, match_score, match_strand, match_phase, match_feat])
                last_lines.insert(0, match_line)
                print("\n".join(last_lines))
#Initialize
                count = 0
                last_chr = this_chr
                last_end = int(this_end)
                last_pos = [int(this_start), int(this_end)]
                last_strand = this_strand
                this_feat = "ID={}-exon-{};Name={};Parent={}".format(feat_name, str(count), feat_name, feat_name)
                line = "\t".join([this_chr, this_type, source, this_start, this_end, this_score, this_strand, this_phase, this_feat])
                last_lines = [line]
                last_name = feat_name
        else:
            match_chr = last_chr
            match_type = "match"
            match_start = str(min(last_pos))
            match_end = str(max(last_pos))
            match_score = '0'
            match_strand = last_strand
            match_phase = '.'
            match_feat = "ID={};Name={}".format(last_name, last_name)
            match_line = "\t".join([match_chr, match_type, source, match_start, match_end, match_score, match_strand, match_phase, match_feat])
            last_lines.insert(0, match_line)
            print("\n".join(last_lines))


if __name__ == "__main__":
    main()
