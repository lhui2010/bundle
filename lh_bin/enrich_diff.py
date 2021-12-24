#!/usr/bin/env python

#TODO: add argparse
#Sample usage:
#python ../enrich_diff.py tmp.441 tmp.445 tmp.diff
#Sample input:
#The following results was generated from bedtools pairtobed. 
#This script would compare those two files by column 11 
#like GNU diff -y tmp.441 tmp.445 -column 11 except that GNU diff do not support diff by column
#==> tmp.441 <==
#M441_chr01      566075  613003  M445_chr01      785040  831972  SYNAL2  M441_chr01      566202  566375  TE_00002369___DNA/Helitron      0       -
#M441_chr01      566075  613003  M445_chr01      785040  831972  SYNAL2  M441_chr01      566441  566557  TE_00012493___DNA/Helitron      0       +
#M441_chr01      566075  613003  M445_chr01      785040  831972  SYNAL2  M441_chr01      566545  566645  TE_00012493___DNA/Helitron      0       +
#M441_chr01      566075  613003  M445_chr01      785040  831972  SYNAL2  M441_chr01      566639  566821  TE_00009470___DNA/DTT   0       +
#M441_chr01      566075  613003  M445_chr01      785040  831972  SYNAL2  M441_chr01      566883  567272  TE_00006172___MITE/DTH  0       -
#M441_chr01      566075  613003  M445_chr01      785040  831972  SYNAL2  M441_chr01      568926  569048  TE_00002124___DNA/Helitron      0       -
#M441_chr01      566075  613003  M445_chr01      785040  831972  SYNAL2  M441_chr01      569141  569378  TE_00002124___DNA/Helitron      0       -
#M441_chr01      566075  613003  M445_chr01      785040  831972  SYNAL2  M441_chr01      570031  570234  TE_00004691___MITE/DTT  0       +
#M441_chr01      566075  613003  M445_chr01      785040  831972  SYNAL2  M441_chr01      570031  570234  570032..570234#MITE/DTT___MITE/DTT      0       .
#M441_chr01      566075  613003  M445_chr01      785040  831972  SYNAL2  M441_chr01      575366  575763  TE_00000883___DNA/DTM   0       -
#
#==> tmp.445 <==
#M441_chr01      566075  613003  M445_chr01      785040  831972  SYNAL2  M445_chr01      784759  785093  TE_00002369___DNA/Helitron      0       +
#M441_chr01      566075  613003  M445_chr01      785040  831972  SYNAL2  M445_chr01      785167  785340  TE_00002369___DNA/Helitron      0       -
#M441_chr01      566075  613003  M445_chr01      785040  831972  SYNAL2  M445_chr01      785406  785522  TE_00012493___DNA/Helitron      0       +
#M441_chr01      566075  613003  M445_chr01      785040  831972  SYNAL2  M445_chr01      785604  785786  TE_00009470___DNA/DTT   0       +
#M441_chr01      566075  613003  M445_chr01      785040  831972  SYNAL2  M445_chr01      785848  786237  TE_00006172___MITE/DTH  0       -
#M441_chr01      566075  613003  M445_chr01      785040  831972  SYNAL2  M445_chr01      787891  788013  TE_00002124___DNA/Helitron      0       -
#M441_chr01      566075  613003  M445_chr01      785040  831972  SYNAL2  M445_chr01      788106  788343  TE_00002124___DNA/Helitron      0       -
#M441_chr01      566075  613003  M445_chr01      785040  831972  SYNAL2  M445_chr01      788996  789199  TE_00004691___MITE/DTT  0       +
#M441_chr01      566075  613003  M445_chr01      785040  831972  SYNAL2  M445_chr01      788996  789199  788997..789199#MITE/DTT___MITE/DTT      0       .
#M441_chr01      566075  613003  M445_chr01      785040  831972  SYNAL2  M445_chr01      792184  792453  TE_00001574___DNA/DTT   0       +

import sys
import re 
import subprocess
import os

#zero based column number
column = 10
oscmd = 'for i in ' + sys.argv[1] + ' ' + sys.argv[2] + '; do awk \'{print $11}\' ${i} | sed \'s/.*#//\' > ${i}.id; done'
os.system(oscmd)
#result = subprocess.run(['awk', '{print $11}', 'tmp.441', '>tmp.441.id'], stdout=subprocess.PIPE)
#result = subprocess.run(['awk', '{print $11}', 'tmp.445', '>tmp.445.id'], stdout=subprocess.PIPE)
result = subprocess.run(['diff', '-y', sys.argv[1] + '.id', sys.argv[2] + '.id'], stdout=subprocess.PIPE)
sysargv3 = str(result.stdout.decode('utf-8')).splitlines()

#    diff -y "tmp.441.id" "tmp.445.id" >"tmp.diff"

#TODO: diff by column
#diff's left column
#8-13
linesA = []
linesB = []
syntag = ''
with open(sys.argv[1]) as fh:
    linesA = fh.readlines() 
linesA = ["\t".join(e.split()[7:13]) for e in linesA]

with open(sys.argv[2]) as fh:
    linesB = fh.readlines() 
syntag = linesB[-1].split()[6]
linesB = ["\t".join(e.split()[7:13]) for e in linesB]


#print(linesB)
#print(len(linesB))
#print(linesB[44])
#print(linesB[45])
#print(linesB[46])
#exit()

linesA_iterator = 0
linesB_iterator = 0
with open(sys.argv[2]) as fh:
    for line in sysargv3:
#Whether to move to next iterator of linesA or liensB
        countA = 1
        countB = 1
#The 
        enriched_A = '\t'.join([''] * 6)
        diff_status = '='
        enriched_B  = '\t'.join([''] * 6)
        mylist = re.split(r'\t+', line.rstrip())
        #print(mylist)
        if(len(mylist) >=3):
            diff_status = mylist.pop(1).strip()
            if diff_status == '>':
                countA = 0
        if line.rstrip().endswith("<"):
            diff_status = '<'
            countB = 0
        if(countA):
        #    print(linesA_iterator)
            enriched_A = linesA[linesA_iterator].rstrip()
            linesA_iterator += 1
        if(countB):
        #    print(linesB_iterator)
            enriched_B = linesB[linesB_iterator].rstrip()
            linesB_iterator += 1
        print("\t".join([syntag, enriched_A, diff_status, enriched_B]))
