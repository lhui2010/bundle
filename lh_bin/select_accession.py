import sys

start_loci = 0
end_loci = 407794900

#Transform MSU7 to RAP2
offset = -3

samp_dict = {}
with open('Sample_Name.with_group.txt') as fh:
#with open('test.header') as fh:
    for line in fh:
        myl = line.rstrip().split('\t')
        samp_dict[myl[0]] = myl[1]

#The list used to store index of favored samples, 1-based
select_list = []
with open(sys.argv[1]) as fh:
    line = fh.readline()
    header_list = line.rstrip().split(',')
    new_header = ''
#get the favored sample colnames as well as corresponding header
    for i in range(1, len(header_list)):
        if header_list[i] in samp_dict:
            select_list.append(i)
            new_header+="\t"+header_list[i]
    print(new_header)
#    print(select_list)
#    exit()
    for line in fh:
        mylist = line.rstrip().split(',')
        loci = int(mylist[0][4:]) + offset
#        mylist.pop(0)
#        mylist.insert(0, str(loci))
#Get the selected genotypes
        is_snp = 1
        if(loci >= start_loci and loci <= end_loci):
            new_line = str(loci)
            for i in select_list:
#Only want SNP now
                if(len(mylist[i]) >1):
                    is_snp = 0
                new_line += "\t" + mylist[i]
#        line = '\t'.join(mylist)
            if is_snp == 1:
                print(new_line)

