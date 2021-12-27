#1												>	1	2466645	2467540	B73_Zm00001d027306_T001	.	-
#1	1	8323799	8332094	B73_Zm00001d027307_T001	0	+	A188G02317-t1	A188G02317-t1	1-to-1	Blast+Nonsyntenic	.	=	1	2523140	2531089	B73_Zm00001d027307_T001	.	+
#1	1	8378959	8384175	B73_Zm00001d049399_T001	0	+	A188G02318-t1	A188G02318-t1	m-to-m.A188	Nonsyntenic	DispersedDuplicates	<						
#1	1	8384908	8388103	B73_Zm00001d027308_T001	0	-	A188G02319-t1	A188G02319-t1	1-to-1	OrthoFinder+Syntenic	.	=	1	2562571	2565147	B73_Zm00001d027308_T001	.	-
#1	1	8429248	8435965	B73_Zm00001d027309_T003	0	-	A188G02320-t1	A188G02320-t1	1-to-1	OrthoFinder+Syntenic	.	=	1	2629289	2635834	B73_Zm00001d027309_T003	.	-
#1												>	1	2650836	2653089	B73_Zm00001d027310_T001	.	-
#1	1	8468222	8471062	B73_Zm00001d027311_T001	0	+	A188G02321-t1	A188G02321-t1	1-to-1	OrthoFinder+Syntenic	.	=	1	2746424	2749031	B73_Zm00001d027311_T001	.	+
#1	1	8510130	8510498	B73_Zm00001d030274_T001	0	-	A188G02322-t1	A188G02322-t1	m-to-m.A188	Nonsyntenic	DispersedDuplicates	<						
#1	1	8511220	8517868	B73_Zm00001d027312_T002	0	+	A188G02323-t1	A188G02323-t1	1-to-1	OrthoFinder+Syntenic	.	=	1	2774054	2780630	B73_Zm00001d027312_T002	.	+
#1	1	8518248	8519587	B73_Zm00001d027313_T001	0	-	A188G02324-t1	A188G02324-t1	1-to-1	OrthoFinder+Syntenic	.	=	1	2781142	2782403	B73_Zm00001d027313_T001	.	-
import sys

start = 0
end = 0

A188_gene_on_B73_loci = ""
B73_gene_on_A188_loci = ""

#The input is actually a kind of enflated bedpe format
class Bedpe():
    """ format the following line into bed information of two sides """
#1	1	8468222	8471062	B73_Zm00001d027311_T001	0	+	A188G02321-t1	A188G02321-t1	1-to-1	OrthoFinder+Syntenic	.	=	1	2746424	2749031	B73_Zm00001d027311_T001	.	+
    def __init__(self, line):
        """Initialize the values"""
        mylist = line.rstrip('\n').split('\t')
        #print(mylist)
        self.gene_left = mylist[7]
        (self.chr_left, self.start_left, self.end_left) = mylist[1:4]
        self.gene_right = mylist[16]
        (self.chr_right, self.start_right, self.end_right) = mylist[13:16]

    def get_right_gene(self):
        return self.gene_right
    
    def get_right_chr(self):
        return self.chr_right

    def get_right_end(self):
        return self.end_right

    def get_right_start(self):
        return self.start_right

    def get_left_gene(self):
        return self.gene_left

    def get_left_chr(self):
        return self.chr_left

    def get_left_start(self):
        return self.start_left

    def get_left_end(self):
        return self.end_left

previous_end_A188 = ""
previous_end_B73 = ""
unsyntenic_gene_B73 = ""
unsyntenic_gene_A188 = ""

with open(sys.argv[1]) as fh:
    for line in fh:
#THis version will not consider gaps
        if 'N{' in line:
            continue
        bedpe = Bedpe(line)
        current_chr_A188 = bedpe.get_left_chr()
        current_start_A188 = bedpe.get_left_start()
        current_end_A188 = bedpe.get_left_end()
        current_gene_A188 = bedpe.get_left_gene()
        current_chr_B73 = bedpe.get_right_chr()
        current_start_B73 = bedpe.get_right_start()
        current_end_B73 = bedpe.get_right_end()
        current_gene_B73 = bedpe.get_right_gene()
        #print("Current gene: ", end ="")
        #print("\t".join([current_gene_A188, current_gene_B73]))
        if '=' in line:
#Use this class to store those variables
#So this is within a syntenic block
            if previous_end_A188 is not "" and previous_end_B73 is not "":
                lociB73 = "\t".join([current_chr_B73, previous_end_B73, current_start_B73])
                lociA188 = "\t".join([current_chr_A188, previous_end_A188, current_start_A188])
                if unsyntenic_gene_A188 is not "":
                    print('{}\t{}'.format(unsyntenic_gene_A188, lociB73))
#Clear geneA188 so that we won't make an error when reading continuous syntenic regions
                if unsyntenic_gene_B73 is not "":
                    print('{}\t{}'.format(unsyntenic_gene_B73, lociA188))
#The same for B73
            unsyntenic_gene_A188 = ""
            unsyntenic_gene_B73 = ""
#So this is within a syntenic block
#            if previous_end_B73 is not empty:
#                print geneA188 and lociB73
#                clear geneA188
            previous_end_A188 = current_end_A188
            previous_end_B73 = current_end_B73
        elif '>' in line:
            unsyntenic_gene_B73 += current_gene_B73 + ","
#            print (line + current_gene_B73)
#            print(line + current_gene_B75)
        elif '<' in line:
            unsyntenic_gene_A188 += current_gene_A188 + ","
#            print (line + current_gene_A188 + "> in line" )
        elif '|' in line:
#            print (line + current_gene_A188 + current_gene_B73)
            unsyntenic_gene_A188 += current_gene_A188 + ","
            unsyntenic_gene_B73 += current_gene_B73 + ","

#        print (line + current_gene_A188 + "ceshi in line" )
