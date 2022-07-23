python get_gene_age.py get_gene_age . all_50.txt.tac2  > gene_pav.txt
selectItem.pl -h nod_genes.txt gene_pav.txt > gene_pav.nodulation.txt
Rscript calc_mrca.R  >gene_pav.nodulation.txt.mrca
sed 's/^\[1\] "//; s/"//g' gene_pav.nodulation.txt.mrca > gene_pav.nodulation.txt.mrca.format
sort -k2,2V gene_pav.nodulation.txt.mrca.format > gene_pav.nodulation.txt.mrca.format.sort
selectItem.pl -k nod_genes.txt.addinfo gene_pav.nodulation.txt.mrca.format.sort  > gene_pav.nodulation.txt.mrca.format.sort.info.txt
