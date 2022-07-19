# From Yan Rui
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import  pandas as pd
import argparse,os
f=[ list(str(i.seq)) for i in SeqIO.parse('/Users/yanrui/Desktop/syb.fas','fasta')]
fs={}
fss=[]
g=[ list(str(i.seq)) for i in SeqIO.parse('/Users/yanrui/Desktop/sar.fas','fasta')]
gs={}
gss=[]


def get_dict(a,li,p):
    lenth = len(a[0])
    for l in range(lenth):
        dic={}
        new={}
        for i in range(len(a)):
            if str(a[i][l]) not in dic:
                dic[str(a[i][l])] = 1
            else:
                dic[str(a[i][l])] += 1
        for i in dic:#取80%以上相似度
            if round(dic[i]/len(a),2) > 0.7 and str(i) != '-':
                # dic[i]=round(dic[i]/len(a),2)
                new[i]= round(dic[i]/len(a),2)
                li[l]=new
                p.append(l)
        # if len(new) != 0 :
        #     li.append([l,new])
        # else:
        #     continue
get_dict(f,fs,fss)
get_dict(g,gs,gss)

diff=sorted(set(fss) & set(gss) )
print('位点','None','postive_select')
with open('/users/yanrui/desktop/3.txt','w') as f:
    for i in diff:
        if gs[i] != fs[i]:
            f.write(str(i)+','+str(gs[i]) + ',' + str(fs[i]) + '\n')
            print(i,gs[i],fs[i])

