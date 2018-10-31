#!/usr/bin/env python
# coding: utf-8
import shlex, subprocess
import re
import sys
from os.path import basename
import os
from Bio import SeqIO

from subprocess import Popen
from itertools import islice

#最多线程数目
max_workers = 50
#最多生成的文件数目
max_partitions= 200

if(len(sys.argv) < 3):
    CRED = '\033[91m'
    CEND = '\033[0m'
    print(CRED + "Usage" + CEND)
    print ("\tpython " + sys.argv[0] + " ref qry ")
    print (CRED+"Example: will split A188 into 500 bp non-overlapping window and align it to B73 using bwa"+CEND)
    print ("\tpython group2fasta.py B73.fasta A188.fasta")
    exit()

qry = sys.argv.pop()
ref = sys.argv.pop()
step = 500

qryseq_dict = SeqIO.to_dict(SeqIO.parse(qry, "fasta"))
refseq_dict = SeqIO.to_dict(SeqIO.parse(ref, "fasta"))

output_file = qry + ref + ".sam"

commands = []
first_prize = 0
for this_qry in qryseq_dict.keys():
    iterr = step
    while iterr <= len(qryseq_dict[this_qry]):
        towrite = qryseq_dict[this_qry][iterr-step:iterr]
        towrite.id = this_qry + "____" +  str(iterr)
        SeqIO.write(towrite, towrite.id , "fasta")
        #$qs 'bwa index Zea_mays.AGPv4.dna_sm.toplevel.fa 
        #bwa mem -t 56 -w 500 -M Zea_mays.AGPv4.dna_sm.toplevel.fa  test.fa '
        #bashCommand = "bwa mem -w 500  -M " + ref + " " + towrite.id + "; rm " + towrite.id 
        if(not first_prize):
            bashCommand = "bwa mem -w 500 -M " + ref + " " + towrite.id  + ">" + output_file + "; rm " + towrite.id
            processes = subprocess.Popen(bashCommand, shell=True)
#            qsub_process = subprocess.Popen(shlex.split(bashCommand), stdout=subprocess.PIPE)
            first_prize = 1
        else:
            bashCommand = "bwa mem -w 500 -M " + ref + " " + towrite.id  + "|grep -v \"^@\" >>" + output_file + "; rm " + towrite.id
        commands.append(bashCommand)
        iterr += step
        #生成1000个文件就开跑bwa
        if(len(commands) >max_partitions or iterr > len(qryseq_dict[this_qry])):
            processes = (Popen(cmd, shell=True) for cmd in commands)
            running_processes = list(islice(processes, max_workers))  # start new processes
            while running_processes:
                for i, process in enumerate(running_processes):
                    if process.poll() is not None:  # the process has finished
                        running_processes[i] = next(processes, None)  # start new process
                        if running_processes[i] is None: # no new processes
                            del running_processes[i]
                            break
            commands = []
       
