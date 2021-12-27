import sys
import re
import logging

#genus dict
genus_dict = "genus2family2order.txt"
genus = {}
with open(genus_dict) as fh:
    for line in fh:
        if(line.startswith("#")):
            continue
        mylist = line.rstrip().split()
        if(len(mylist) < 3):
#            logging.warning(mylist)
            mylist.append('')
#        genus[mylist[0]] = "\t".join([mylist[2], mylist[1], mylist[0]])
        genus[mylist[0]] = line.strip('\n')

#:%s/<a href="//g
#:%s/" class="vis-tooltip">/\t/g
#:%s/<.\{-\}>/\t/g

pattern1 = r'<a href="'
pattern2 = r'" class="vis-tooltip">'
pattern3 = r'<.*?>'


header = "\t".join(["Order", "Family", "Genus", "Species name", "Common name", "Genome size", "Author", "Journal", "Title", "Doi"])
with open(sys.argv[1]) as fh:
    for line in fh:
        line = re.sub(pattern1, '\n', line)
        line = re.sub(pattern2, '\t', line)
        line = re.sub(pattern3, '\t', line)
        line = re.sub(r'\t+', '\t', line)
#        if(line.strip() != ""):
#            print(line)
#            exit()
        lines = line.split('\n')
        for l in lines:
            if l.strip() == "":
                continue
            mylist = l.split('\t')
            (doi, sp_name, sp_name2, common_name, size, author, title, journal) = mylist[:8]
#0                                          1                           2               3                       4                   5                       6   ..7
#http://dx.doi.org/10.1111/tpj.14744    Populus ilicifolia    Populus ilicifolia     (African poplar)    genome size: 444 Mbp    Chen Z et al. (2020)    Survival in the Tropics despite isolation, inbreeding and asexual reproduction: insights from the genome of the world's southernmost poplar (Populus ilicifolia).    Plant J. 2020 Mar 13, Epub ahead of print
            logging.warning(sp_name.split())
            (genus_name, species_name) = sp_name.split( maxsplit=1)

            if('20' in journal):
                (journal, pub_date) = journal.split('20', maxsplit = 1)
                pub_date = "20" + pub_date
            else:
                (author, pub_date) = author.split('(', maxsplit = 1)
                pub_date = pub_date.rstrip(')')

            if genus_name in genus:
                genus_info = genus[genus_name]
            else:
                genus_info = ["\t"]*2

            common_name = common_name.lstrip()
            common_name = common_name.lstrip('(')
            common_name = common_name.rstrip(')')
            new_line = "\t".join([sp_name, common_name, journal, pub_date, author, title, size,  doi, genus_info])
            print(new_line)


