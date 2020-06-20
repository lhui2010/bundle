
import sys
#Hic.fastq.gz.counts_GATC.20g1	1	79244	1	W	tig00005006|arrow_np1212	1	79244	+
#Hic.fastq.gz.counts_GATC.20g1	79245	79344	2	U	100	contig	yes	map
#Hic.fastq.gz.counts_GATC.20g1	79345	123128	3	W	tig00005007|arrow_np1212	1	43784	-
#Hic.fastq.gz.counts_GATC.20g1	123129	123228	4	U	100	contig	yes	map
#Hic.fastq.gz.counts_GATC.20g1	123229	195902	5	W	tig00000193|arrow_np1212	1	72674	-
#Hic.fastq.gz.counts_GATC.20g1	195903	196002	6	U	100	contig	yes	map
#Hic.fastq.gz.counts_GATC.20g1	196003	258290	7	W	tig00000184|arrow_np1212	1	62288	-
#Hic.fastq.gz.counts_GATC.20g1	258291	258390	8	U	100	contig	yes	map
#Hic.fastq.gz.counts_GATC.20g1	258391	345821	9	W	tig00000180|arrow_np1212	1	87431	+
#Hic.fastq.gz.counts_GATC.20g1	345822	345921	10	U	100	contig	yes	map

contig_list = []
order_list = []
order_index = -1
contig_index = 0
last_tag = ""

with open(sys.argv[1]) as fh:
    for line in fh:
        mylist = line.rstrip().split()
        if (mylist[4] == "U"):
            continue
        else:
            contig_index +=1
        contig_list.append(mylist[5])
        order_abbrev = int(mylist[8]) * contig_index
        if(last_tag != mylist[0]):
            order_index += 1
            order_list[order_index] = ""
        if(last_tag !=
        
        order_list[order_list] += formated_order

        
        last_tag = mylist[0]
