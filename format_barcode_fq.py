#!/usr/bin/env python
import argparse
import re
import gzip
#import sys

#Args
parser = argparse.ArgumentParser(description='Trans form barcoded.fq.gz from Longranger to \
fastq needed by SLR-superscaffolder')
parser.add_argument('INPUT', type=str, nargs = 1,
                    help='input fastq file')
parser.add_argument('--version', action='version', version='%(prog)s 0.1')
args = parser.parse_args()

inputs = args.INPUT[0]

"""
Input fastq: 
@E00591:54:H53WFCCXY:4:1105:16995:72983
AAAGATTTCTATGAAAGCTTTCTGTTTTGAGCGCTATACATTTTTTAATAAGGGAGCAAATAACAGAAAAAACAAACATACACACACAAATGAGGGATTATGAAATCACTTACTCCCTGTTGTGCCA
+
FJ-7AJ-FFA-<FJJJFF<FJ-<F--FAJ7F-F<-FF--F<FF-FA<7-<<-<-<A7777--7<F--AF-7FFFJJJF-7F<AA-AFFJJA7<A----7777FJJA-FJJ<FA---A-J--7A--A7
@E00591:54:H53WFCCXY:4:1105:16995:72983
AGTGTTCTGCCGGGTAGAAGGAGGACAGATAGCTGCAATGATCTTCCTGTTTGTCACCATGATAGTTTATCTCATTAGTGCGTTGGTTTGCGTAAAGTTATGGAGGCATGAGGCAGCTCGGAGACATTGAGAGTATATGGAACAACAGGA
+
A<A-AFJA-F<A-7-7FF-FJ<77--<<7FFF<<--7-FJAF-AA<7F-<A<F7FJJJJFFF<-FAAJ<FF--7A<<F<J7-----<-FA<-A7FJ7A-A-77F<<-AJAF---A-A<)-7-))--A-77<7-<---7AF-77---7A7-

Output fastq_1.fq:
@ST-E00575:76:H3KLKCCXY:7:2111:14093:13246#1113133143412121/1
ATATCGCTGGTGTTTAGATGGGCGGTACCCTTGGAAGTTGAGCAGGAATTTATCCCGGAGGCTTCTCCAGGAGGCAATGGACAGCGGAGGCAACCTGGTGAACCAGGTGAGAGCTGGGCCGTGGAGG
+
JFJJ-7<-7A7AA-7-<-<--A<FA--F7FA-77-FF7A-<-7F<JJ-<AJFA--7A--7<JA7--A<FJJAJ7A<AAAA-A<-7AAF7F<FJ7<--7FA7<7-A-)<A7A-<)-<<)-AA7)<-7A
Output fastq_2.fq:
@ST-E00575:76:H3KLKCCXY:7:2111:14093:13246#1113133143412121/2
CTGGAGGCCAAGCGCCAGCGCGTGTCCGCACTAGCCAAGGTGCGGCAAATGATACGCGACAAAGAGCAGAAGGCCCGGGACCTCGAGCGGGAGATCGCGCTGATGCAGCGCGAAGGCCACCTCGGCCCGCCGCATGGGCCACCCCTGCAA
+
AAFFFFFJ<FFF7FF<<JF<<<-<--FAFF<JF<FJ-7F77<FJ<-<AAJFA<FFJF-F-7AJ<FA-77A-<AAJF-7<FJAA7<-7AAF-AF--7-FF7J<F7FF7-7-A-)A))7FJ-A)FA<)-))7)-7))<-7<-7<F7<A<))-
"""

ATCGN_dict = {'A':'1', 
    'T':'2',
    'C':'3',
    'G':'4',
    'N':'0'}
def digit_encode(sequence):
    res = ''
    for i in sequence:
        res += ATCGN_dict[i]
    return res

barcode_list = []
count = 0
if(re.search(r'gz$', inputs)):
    fh = gzip.open(inputs)
else:
    fh = open(inputs)

#with eval(open_string) as fh, \
with open (inputs + "_1.fq", "w") as fhpe1, \
    open (inputs + "_2.fq", "w") as fhpe2:
    line = fh.readline().decode('ASCII')
    while line:
        """Read and setting flag"""
        header = line
        line_rest = ''
        for i in range(3):
            line_rest += fh.readline().decode('ASCII')
        """TestWhether this reads is barcoded"""
        if 'BX:' in header:
            """Change seperater and add suffix"""
            #header = re.sub(r'\sBX.*:', '#', header)
            header = re.sub(r'-\d', '', header).rstrip()
            (left, right) = re.split(r'\sBX.*:',  header)
            right = digit_encode(right)
            header = left + "#" + right
            #print(left)
            #print(right)
            #print(new_right)
            #exit()
            if(header in barcode_list):
                header += "/2\n"
                fhpe2.write(header + line_rest)
                count += 1
                barcode_list = []
            else:
                barcode_list.append(header)
                header += "/1\n"
                fhpe1.write(header + line_rest)
        line = fh.readline().decode('ASCII')
fh.close()
