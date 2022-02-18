#!/usr/bin/env python

# -*- coding = utf-8 -*-
# @Time : 2022/2/18 9:37
# @Author : 严慧
# @File : count_GC%.py
# @Software : PyCharm

import sys
from Bio import SeqIO
from Bio.SeqUtils import Seq
from Bio.SeqUtils import GC
#GC 计算G+C含量，返回百分比
from Bio.SeqUtils import GC123
#计算GC含量，返回4个百分比：合计和三个密码子位置的GC含量
import textwrap
import argparse





def count_GC(seq):
    GC1 = []
    GC2 = []
    GC3 = []
    for i in range(1,len(seq)+1):
        if i%3 == 1:
            GC1.append(str(seq[i-1]))
        elif i%3 == 2:
            GC2.append(str(seq[i-1]))
        else:
            GC3.append(str(seq[i-1]))
    GC1_seq = Seq("".join(GC1))
    GC2_seq = Seq("".join(GC2))
    GC3_seq = Seq("".join(GC3))
    return [round(GC(GC1_seq), 2), round(GC(GC2_seq), 2), round(GC(GC3_seq), 2), round(GC(seq), 2)]


def main():
######################
### Argument parse ###
######################
    prog_name = "count_GC"
    usage = r'''
计算一段序列第三位密码子的GC含量
比如ACGTCGAAATTT，算出来的GC%就是50%，如果是ACGTCGACCTCC，GC%就是100%
输入fasta文件，输出一个table，左边是每段fasta名称，右边是fasta的三位密码子的GC%
'''
#use of textwrap allowed custome format like "\t\tusage"
    parser = argparse.ArgumentParser(
        prog=prog_name, 
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(usage), 
        epilog="")
    parser.add_argument("fasta", help="input fasta file")
    parser.add_argument("gc_table", help="output gc table file")
    args = parser.parse_args()  
    in_fasta = args.fasta
    out_table = args.gc_table


    with open(out_table,"w") as fw:
        fw.write("%s\t%s\t%s\t%s\t%s\n"%("gene_name","GC1","GC2","GC3","GC_total"))
        buff = ''
        for rec in SeqIO.parse(in_fasta, "fasta"):
            # Bio.SeqIO.parse(文件句柄或字符串形式的文件名，文件格式的小写字符串)，将序列文件转换为SeqRecords的迭代器
            GC_result = GC123(rec.seq)
            GC_all = round(GC_result[0], 2)
            GC_1 = round(GC_result[1], 2)
            GC_2 = round(GC_result[2], 2)
            GC_3 = round(GC_result[3], 2)
            # GC_num = count_GC(rec.seq)
            # GC_1 = GC_num[0]
            # GC_2 = GC_num[1]
            # GC_3 = GC_num[2]
            # GC_all = GC_num[3]
            buff += ("%s\t%s\t%s\t%s\t%s\n" % (rec.id, str(GC_1), str(GC_2), str(GC_3), str(GC_all)))
        fw.write(buff)

    # file = open(inputfile, "r")
    # with open(outputfile, "w") as f:
    #     f.write("name", "\t", "GC%", "\t", "GC%for1", "\t", "GC%for2", "\t", "GC%for3")
    #     for num, value in enumerate(file):
    #         if value[0] == ">":
    #             name = value.replace(">", "")
    #             content = next(file)
    #             content_1 = content[::3]
    #             content_2 = content[1::3]
    #             content_3 = content[2::3]
    #             G = content.count("G")
    #             C = content.count("C")
    #             GC = G + C
    #             GCP = GC / len(content) * 100
    #             G1 = content_1.count("G")
    #             C1 = content_1.count("C")
    #             GC1 = G1 + C1
    #             GCP1 = GC1 / len(content_1) * 100
    #             G2 = content_2.count("G")
    #             C2 = content_2.count("C")
    #             GC2 = G2 + C2
    #             GCP2 = GC2 / len(content_2) * 100
    #             G3 = content_3.count("G")
    #             C3 = content_3.count("C")
    #             GC3 = G3 + C3
    #             GCP3 = GC3 / len(content_3) * 100
    #             f.write(name, "\t", GCP, "\t", GCP1, "\t", GCP2, "\t", GCP3)

if __name__ == '__main__':
    main()
