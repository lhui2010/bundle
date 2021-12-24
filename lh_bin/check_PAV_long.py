#!/usr/bin/env python
import re
import argparse
from Bio import SeqIO
import sys
import subprocess
import multiprocessing
from pprint import pprint


class GenewiseResult:
    """A simple object to store genewise results. Required
    parameters are pretty(str) cds(str) pep(str) gff(str) 
    and pesudo(bool)"""
    def __init__ (self, pretty, cds, pep, gff, pseudo):
        self.pretty = pretty
        self.cds = cds
        self.pep = pep
        self.gff = gff
        self.pseudo = pseudo

#    def format_gff(gff, cut_start, cut_end , strand = '+', qry_name):
#    """Change coordinates, cut_start and cut_end is 1-based seq from big seq, 
#    and strand specifies whether this gff was on an antisense strand.
#    """
#        #Change qry_name

#    def format_cds():
#    """Genewise usually missed TAA at the end. Here we will add TAA at the 
#    end of CDS sequences if TAA exists on genome sequence.
#    """
    #genwise_result.pretty = str
    #genwise_result.cds = str
    #genwise_result.pep = str
    #genwise_result.gff = str
    #genwise_result.pseudo = bool 

def pretty2seq(pretty_str):
    "Transform pretty format to cds format"
    pretty_cds = ''
    pretty_list = pretty_str.split('\n')
    for i in range(10, len(pretty_list), 8):
        pretty_cds += pretty_list[i].strip() + "\n"
    pretty_cds = re.sub(r'[\- ]', '', pretty_cds)
    pretty_cds = pretty_cds.replace('!', '*')
    return pretty_cds

#Args
def checkPAV(bedline, qry_dict, ref_dict, pretty_dict):
    off_set = 0
    split = '_'
#1_RaGOO 5336153 5342648 B73_Zm00001d027233_T001 60      -
    bedlist = bedline.rstrip().split()
    (ref_chr, ref_start, ref_end, qry_name, score, ref_strand) = bedlist[0:6]
#Avoid name discordance
    ref_chr = ref_chr.replace('_RaGOO', '')
    flank_start = int(ref_start) - off_set
    flank_end = int(ref_end)  + off_set
    python_start = flank_start - 1
    python_end = flank_end 
    #ref_seq = ref_dict[ref_chr].seq[(int(ref_start)-1):int(ref_end)]
    ref_seq = ref_dict[ref_chr].seq[python_start:python_end]
    if(ref_strand == '-'):
        ref_seq = ref_seq.reverse_complement()
#Generate Qry.fasta and Ref.fasta
#    SeqIO.write(xx)
#$genewise Zm00001d042922_T002.pep target.fa >Zm00001d042922_T002.fa.vs.target.gene_wise
#    genewise = 'genewise %s'%(qry_name, target_name
    qry_seq = qry_dict[qry_name].seq#.__str__
    qry_len = len(qry_seq)

    while(qry_name in pretty_dict):
        qry_name += "DUP"

    qry_fasta = qry_name + '.fa'
    ref_fasta = ref_chr + ref_start + ref_end + '.fa'
    with open(qry_fasta, 'w') as fqry, \
    open(ref_fasta, 'w') as fref:
    #    SeqIO.write(qry_seq, fqry, "fasta")
    #    SeqIO.write(ref_seq, fref, "fasta")
        fqry.write(">%s\n%s"%(qry_name, qry_seq))
        fref.write(">%s\n%s"%(ref_chr + ref_start + ref_end, ref_seq))
    wise_pretty = subprocess.check_output("genewise -pseudo -pretty '%s' '%s' "%(qry_fasta, ref_fasta) , \
    shell=True).decode("utf-8")
    
#genewise -gff -pseudo -cdna -trans test_qry.fa test_ref.pseudo.fa
    wise_output = subprocess.check_output("genewise -gff -cdna -pseudo -trans '%s' '%s' "\
    %(qry_fasta, ref_fasta) , shell=True)
    wise_output = wise_output.decode("utf-8")
#    rm = 'rm %s %s'(qry_fasta, ref_fasta)
#    os.system(rm)
#    print(wise_output)
    
    (wise_pep, wise_cds, wise_gff, tmp) = wise_output.split(sep='//')

    if 'pseudo' in wise_pep:
        wise_pseudo = True
    else:
        wise_pseudo = False

    predict_len = 0
#    if not wise_pseudo:
    (predict_id, predict_seq) = wise_cds.strip().split('\n', 1)
    predict_cds = predict_seq.replace('\n', '')
    predict_len = len(predict_cds)/3
    wise_pep = re.sub(r'>.*?\n', ">"+qry_name+"\n", wise_pep)
#    else:
#TODO export PRETTY to PEP
    if wise_pseudo:
        wise_pep = ">" + qry_name + "\n" + pretty2seq(wise_pretty)
#Mark pseudo gene
    if wise_pseudo or '*' in wise_pep:
        qry_name += '.0'

    wise_pep = re.sub(r'>.*?\n', ">"+qry_name+"\n", wise_pep)
    wise_cds = re.sub(r'>.*?\n', ">"+qry_name+"\n", wise_cds)

    wise_pseudo = "\t".join([str(wise_pseudo), "%.2f"%(predict_len/qry_len), str(predict_len), str(qry_len)])

## Foramt GFF
###Format stopcodon
    gff_lines = wise_gff.split('\n')
    if(gff_lines[-1].strip() == ''):
        gff_lines.pop()
    if(gff_lines[0].strip() == ''):
        gff_lines.pop(0)
    mylocilist = []
    #print(gff_lines[-1])
    myl = gff_lines[-1].split()
    predict_end = int(myl[4])
    if(ref_seq[predict_end:predict_end+3] == 'TAG' or 
        ref_seq[predict_end:predict_end+3] == 'TAA' or
        ref_seq[predict_end:predict_end+3] == 'TGA'):
        extend_stop_codon = True
    myl[4] = str(predict_end + 3)
    gff_lines[-1] = "\t".join(myl)
    wise_cds = wise_cds.strip() + str(ref_seq[predict_end:predict_end+3]) +"\n"
    #print(gff_lines[-1])

### Format Tag
    for i,g in enumerate(gff_lines):
        myl = g.split()
        myl[0] = ref_chr
        #myl[-1] = "ID=" + qry_name
        if(myl[2] == 'match'):
            myl[2] = 'mRNA'
            myl[-1] = "ID=" + qry_name + ";Name="+qry_name
        elif(myl[2] == 'intron'):
            myl[2] = 'INTRON'
            myl[-1] = "ID=" + qry_name + ".INTRON;Parent="+qry_name
        elif(myl[2] == 'cds'):
            myl[2] = 'CDS'
            myl[-1] = "ID=" + qry_name + ".CDS;Parent="+qry_name
### Format coordinate
        if(ref_strand == '+'):
            myl[3] = str(int(ref_start) - 1 + int(myl[3]))
            myl[4] = str(int(ref_start) - 1 + int(myl[4]))
        else:
#### Reverse
            myl[3] = str(int(ref_end) - int(myl[3]) + 1)
            myl[4] = str(int(ref_end) - int(myl[4]) + 1)
            (myl[3], myl[4]) = (myl[4], myl[3])
            myl[6] = '-'
        gff_lines[i] = "\t".join(myl)
        #print(ref_start)
        #print(myl[3]+"\t"+myl[4])
        #print(gff_lines[i])
#### Transform coordinate
    #print(gff_lines[-1])
    wise_gff = "\n".join(gff_lines) + "\n"

    pretty_dict[qry_name] = GenewiseResult(wise_pretty, wise_cds, wise_pep, wise_gff, wise_pseudo)

#    wise_cds = '>' + wise_cds


if __name__ == '__main__':
#Program
    parser = argparse.ArgumentParser(description='''Regex based exact sequence alignment tool\nResult is 1-based
    Eg: check_PAV.py A188.B73.ortho_B73_unique.cds.bed A188.B73.ortho_B73_unique.pep A188.genome
    ''')
    parser.add_argument('BED', type=str, nargs = 1,
                        help='bed file of candidant PAV regions of query pep on reference genome')
    parser.add_argument('PEP', type=str, nargs = 1,
                        help='fasta of query peptide sequences')
    parser.add_argument('REF', type=str, nargs = 1,
                        help='fasta of reference sequences to be checked with PAV')
    parser.add_argument('-t', '--threads', type=str, nargs = '?',
                        help='fasta of reference sequences to be checked with PAV')
    parser.add_argument('--version', action='version', version='%(prog)s 0.1')
    args = parser.parse_args()

#Prepare fasta
    sys.stderr.write("Reading fasta sequences..\n")
    qry_dict = SeqIO.to_dict(SeqIO.parse(args.PEP[0], "fasta"))
    ref_dict = SeqIO.to_dict(SeqIO.parse(args.REF[0], "fasta"))

    with open(args.BED[0]) as fbed:
        bed_list =fbed.readlines()
    result_dict = {}
    for i in range(0, len(bed_list)):
        checkPAV(bed_list[i], qry_dict, ref_dict, result_dict)
    #pool = multiprocessing.Pool(int(args.threads))
    #pool.map(checkPAV, bed_list)
    #pool.close()
    #pool.join()
    prefix = args.BED[0]
    fname_pretty = prefix + ".pretty"
    fname_cds = prefix + ".cds"
    fname_gff = prefix + ".gff"
    fname_pep = prefix + ".pep"
    fname_pseudo = prefix + ".pseudo"
    with open(fname_pretty, 'w') as fprettty, \
    open(fname_cds, 'w') as fcds, \
    open(fname_gff, 'w') as fgff, \
    open(fname_pep, 'w') as fpep, \
    open(fname_pseudo, 'w') as fpseudo:
        for k in result_dict:
            fprettty.write(result_dict[k].pretty)
            fcds.write(result_dict[k].cds)
            fgff.write(result_dict[k].gff)
            fpep.write(result_dict[k].pep)
            fpseudo.write(k + "\t" + str(result_dict[k].pseudo) + "\n")
