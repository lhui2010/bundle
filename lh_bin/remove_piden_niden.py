#!/usr/bin/env python

import sys
import os

# The blast command and outfmt
# blastp -db Andira_inermis_Pap.pep -query Andira_inermis_Pap.pep -out Andira_inermis_Pap.pep.yangya.bln -evalue 1e10 -outfmt "6 qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore" -num_threads 10
#remove larger gene families
rawblast = sys.argv[1]
filteredblast = rawblast + '.filter'
if not os.path.exists(filteredblast):
    block = [] # store a list of query,hit pairs with the same query
    infile = open(rawblast,"r")
    outfile = open(filteredblast,"w")
    last_query,last_hit = "",""
    for line in infile:
        spls = line.split("\t")
        if len(spls) < 3: continue
        query,hit,pident,nident = spls[0],spls[2],float(spls[5]),int(spls[6])
        if query == hit or pident < 20.0 or nident < 50: continue
        if query == last_query or block == []:
            if (query,hit) not in block:
                block.append((query,hit)) # add tuple to block
        else:
            if len(block) < 10:
                for i in block:
                    outfile.write(i[0]+"\t"+i[1]+"\n")
            block = []#reset and initiate the block
            block.append((query,hit))    
        last_query,last_hit = query,hit
    # process the last block
    if len(block) < 10:
        for i in block:
            outfile.write(i[0]+"\t"+i[1]+"\n")
    block = []#reset and initiate the block
    infile.close()
    outfile.close()
