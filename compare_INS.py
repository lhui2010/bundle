#!/usr/bin/env python

#This is a script used to work the output of enrich_diff.py

#SYNAL1	M441_chr01	565790	566066	TE_00002369___DNA/Helitron	0	+	=	M445_chr01	784759	785093	TE_00002369___DNA/Helitron	0	+
#SYNAL2	M441_chr01	566202	566375	TE_00002369___DNA/Helitron	0	-	=	M445_chr01	784759	785093	TE_00002369___DNA/Helitron	0	+
#SYNAL2	M441_chr01	566441	566557	TE_00012493___DNA/Helitron	0	+	|	M445_chr01	785167	785340	TE_00002369___DNA/Helitron	0	-
#SYNAL2	M441_chr01	566545	566645	TE_00012493___DNA/Helitron	0	+	=	M445_chr01	785406	785522	TE_00012493___DNA/Helitron	0	+
#SYNAL2	M441_chr01	566639	566821	TE_00009470___DNA/DTT	0	+	=	M445_chr01	785604	785786	TE_00009470___DNA/DTT	0	+
#SYNAL2	M441_chr01	566883	567272	TE_00006172___MITE/DTH	0	-	=	M445_chr01	785848	786237	TE_00006172___MITE/DTH	0	-
#SYNAL2	M441_chr01	568926	569048	TE_00002124___DNA/Helitron	0	-	=	M445_chr01	787891	788013	TE_00002124___DNA/Helitron	0	-
#SYNAL2	M441_chr01	569141	569378	TE_00002124___DNA/Helitron	0	-	=	M445_chr01	788106	788343	TE_00002124___DNA/Helitron	0	-
#SYNAL2	M441_chr01	570031	570234	TE_00004691___MITE/DTT	0	+	=	M445_chr01	788996	789199	TE_00004691___MITE/DTT	0	+
#SYNAL2	M441_chr01	570031	570234	570032..570234#MITE/DTT___MITE/DTT	0	.	=	M445_chr01	788996	789199	788997..789199#MITE/DTT___MITE/DTT	0	.
import sys

#case flag:
#   0: last line is == and current line is ==
#   1: last line is == and current line is >/<
#   2: last line is </> and current line is ==
ratio_cutoff = 0.5
distance_cutoff = 300
last_type = '='
flag = 0
with open(sys.argv[1]) as fh:
    for line in fh:
        myline = line.rstrip()
        mylist = myline.split('\t')
        current_type = mylist[7]
        if(current_type == '>' or current_type == '<'):
            flag = 1
            continue
        elif(current_type == '='):
            (current_chr, current_start1, current_end1) = mylist[1:4]
            (current_start2, current_end2) = mylist[9:11]
            if(flag and current_chr == last_chr):
                distance1 = abs(int(current_start1) - int(last_end1))
                distance2 = abs(int(current_start2) - int(last_end2))
                if(distance1 == 0 and distance2 == 0):
                    ratio = 1
                else:
                    ratio = min(distance2, distance1) / max(distance1, distance2)
                #print("\t".join([current_chr, current_start1, current_end1, distance1, distance2]))
                if( ratio <= ratio_cutoff and  abs(distance2 - distance1) >= distance_cutoff):
                    print ("*", end = "")
                print("\t".join([current_chr, current_start1, current_end1, str(distance1), str(distance2), str(ratio)]))
            flag = 0
            last_end1 = current_end1
            last_end2 = current_end2
            last_type = current_type
            last_chr = current_chr
            
