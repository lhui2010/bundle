import sys

ortho_dict = {}
with(open(sys.argv[1])) as fh:
    for l in fh:
        if(l.startswith(">")):
            newl = l.rstrip().replace(">", "")
            ortho_gene = newl
            (original_gene, suffix) = newl.split('-R')
            if original_gene not in ortho_dict:
                ortho_dict[original_gene] = ortho_gene
            else:
                ortho_dict[original_gene] += "," + ortho_gene

for k in ortho_dict:
    print(k + "\t" + ortho_dict[k])
