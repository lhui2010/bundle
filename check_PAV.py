#!/usr/bin/env python
import re
import argparse
from Bio import SeqIO
import sys
import subprocess

#Args
def checkPAV(bedline, qry_dict, ref_dict):
    split = '_'
#1_RaGOO 5336153 5342648 B73_Zm00001d027233_T001 60      -
    bedlist = bedline.rstrip().split()
    (ref_chr, ref_start, ref_end, qry_name, score, ref_strand) = bedlist[0:6]
#Avoid name discordance
    ref_chr = ref_chr.replace('_RaGOO', '')
    ref_seq = ref_dict[ref_chr].seq[(int(ref_start)-1):int(ref_end)]
    if(ref_strand == '-'):
        ref_seq = ref_seq.reverse_complement()
#Generate Qry.fasta and Ref.fasta
#    SeqIO.write(xx)
#$genewise Zm00001d042922_T002.pep target.fa >Zm00001d042922_T002.fa.vs.target.gene_wise
#    genewise = 'genewise %s'%(qry_name, target_name
    qry_seq = qry_dict[qry_name].seq#.__str__

    qry_fasta = qry_name + '.fa'
    ref_fasta = ref_chr + ref_start + ref_end + '.fa'
    with open(qry_fasta, 'w') as fqry, \
    open(ref_fasta, 'w') as fref:
    #    SeqIO.write(qry_seq, fqry, "fasta")
    #    SeqIO.write(ref_seq, fref, "fasta")
        fqry.write(">%s\n%s"%(qry_name, qry_seq))
        fref.write(">%s\n%s"%(ref_chr + ref_start + ref_end, ref_seq))
    wise_prediction = subprocess.check_output("genewise -pretty '%s' '%s' "%(qry_fasta, ref_fasta) , shell=True)
#genewise -gff -pseudo -cdna -trans test_qry.fa test_ref.pseudo.fa
    if('!' in wise_output):
        pseudo = True
    wise_output = subprocess.check_output("genewise -gff -cdna -pseudo -trans '%s' '%s' "%(qry_fasta, ref_fasta) , shell=True)
#    seq1 = subprocess.check_output("blastdbcmd -db '%s' -entry '%s' -range '%s'-'%s'" % (genomeFile+spliter+"db", entry, int(p_start), int(p_end)), shell=True)
    wise_output = wise_output.decode("utf-8")
#    rm = 'rm %s %s'(qry_fasta, ref_fasta)
#    os.system(rm)
    #print(wise_output)
    (wise_cds, wise_gff, tmp) = wise_output.split(sep='//')
    print(wise_cds)
    print(wise_gff)

#    fh_cds.write('>' + ortholog_prefix + gene + '\n' + cds_dict[gene].__str__() + '\n')
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
    checkPAV(bed_list[0], qry_dict, ref_dict)
#    pool = multiprocessing.Pool(int(args.threads))
#    pool.map(checkPAV,l)
#    pool.close()
#    pool.join()

    
