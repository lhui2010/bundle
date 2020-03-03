#!/usr/bin/env python

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



ortholog_prefix = ""

def parse_hsp(line):
    this_array = re.split('\(|\)|-', line)
    return this_array

#Zm00001d015434_T001|scaffold35,8459160,:6592798..6593478|+|gene cover:681(100%)|score:336|rank:1
#HSP_ID[1]:(6592798-6593478);query:(1-681); pid: 98.6784

def compile_bed(line, hsp, genome_dict):
    this_array = re.split('\||:|\.\.', line)
    gene_id = re.sub(r'-t\d+$', '', this_array[0])
    contig_name = re.sub(r',.*', '', this_array[1])
    start = this_array[2]
    end = this_array[3]
    strand = this_array[4]
    score = this_array[8]
#    contig_name = this_array[1] + '|' + this_array[2] + '|' + this_array[3] + '|' + this_array[4]
#    start = this_array[5]
#    end = this_array[6]
#    strand = this_array[7]
#    score = this_array[11]
    bed_compile = "\t".join((contig_name, start, end, ortholog_prefix +gene_id, score, strand))
    feat= "ID=" + gene_id + ";Name=" + gene_id +";\n"
    gff_compile = "\t".join((contig_name, ".", "gene", start, end, ".", strand, ".", feat))
    mRNA_id=gene_id + "-t1"
    mRNAfeat= "ID=" + mRNA_id + ";Parent=" + gene_id +";Name=" + mRNA_id + ";\n"
    gff_compile += "\t".join((contig_name, ".", "mRNA", start, end, ".", strand, ".", mRNAfeat))

##Zm00001d009939_T001|001704F|arrow|arrow|pilon:30526..31453|+|gene cover:618(100%)|score:308.975|rank
##HSP_ID[2]:(31217-31453);query:(382-618); pid: 100
#HSP_ID[3]:(30895-31050);query:(237-392); pid: 99.359
#HSP_ID[1]:(30526-30763);query:(1-238); pid: 100
    hsp_string = ""
    CDSfeat = "ID=" + mRNA_id+":cds;Parent=" + mRNA_id + ";\n"
    cds_compile = ""
    for h in hsp:
        (hsp_id, rstart, rend, qry, qstart, qend, pid) = parse_hsp(h)
        h_seq = genome_dict[contig_name].seq[int(rstart)-1:int(rend)]
        if(strand == '-'):
            hsp_string += h_seq.reverse_complement()
            cds_compile += "\t".join((contig_name, ".", "CDS", rstart, rend, ".", strand, ".", CDSfeat)) 
        else:
            cds_compile = "\t".join((contig_name, ".", "CDS", rstart, rend, ".", strand, ".", CDSfeat)) +cds_compile
            hsp_string = h_seq + hsp_string
    gff_compile += cds_compile
    return (gene_id, bed_compile, gff_compile, hsp_string)


#strip suffix of fasta headers like ,123,fragm6c1
def rename_dict(cur_dict):
    return 1


def main():

    parser = argparse.ArgumentParser(description='Extract cds, bed, and gff from genblasta result')
    parser.add_argument('GBLAST', type=str, nargs = 1,
                        help='out.gblast')
    parser.add_argument('FA', type=str, nargs = 1,
                        help='Genome.fasta')
    parser.add_argument('--version', action='version', version='%(prog)s 0.1')


    args = parser.parse_args()

#if args.pos_arg > 10:
#        parser.error("pos_arg cannot be larger than 10")


    fgenblast=args.GBLAST[0]
    ffasta=args.FA[0]

    raw_dict = SeqIO.to_dict(SeqIO.parse(ffasta, "fasta"))            
    genome_dict = {}
#rename genome_dict
    for old_key in raw_dict.keys():
        new_key = re.sub(r',.*', '', old_key)
        genome_dict[new_key] = raw_dict[old_key]

    bed_dict={}
    cds_dict={}
    gff_dict={}

    with open (fgenblast) as fh:
        for line in fh:
            if(re.search(r'rank:1$', line)):
#process bed
                rank_line = line
                line = fh.readline()
                hsp = []
                while(re.match('HSP', line)):
#process cds    
                    hsp.append(line)
                    line = fh.readline()
                (gene_id, bed_compile, gff_compile, cds_compile) = compile_bed(rank_line, hsp, genome_dict)
                bed_dict[gene_id] = bed_compile
#TODO:Create a sequence object here, with name in it ,so it's easier to write to file
                cds_dict[gene_id] = cds_compile
                gff_dict[gene_id] = gff_compile

    with open(fgenblast+".bed", 'w') as fh_bed, \
        open (fgenblast + ".cds", 'w') as fh_cds, \
        open (fgenblast + ".gff", 'w') as fh_gff:
        for gene in bed_dict:
            fh_bed.write('{}\n'.format(bed_dict[gene]))
            fh_gff.write('{}'.format(gff_dict[gene]))
            fh_cds.write('>' + ortholog_prefix + gene + '\n' + cds_dict[gene].__str__() + '\n')

main()
