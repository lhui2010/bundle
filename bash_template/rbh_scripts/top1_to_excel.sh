#rm all_join.txt*
#  for i in Medicago_truncatula.*top1; do awk '{print $2"\t"$1}' $i > $i.ortho; done
#  for i in *.ortho; do perl ../../02.nodulation_genes/dup2col.pl $i >> all_join.txt; done
#  python ../../02.nodulation_genes/stat2matrix.py all_join.txt
sed -i "s/$/\$/" all_join.txt.xls
perl ../../03.blast_rbh/count_matrix.pl all_join.txt.xls > all_join.txt.xls.count
sed -i "s/mRNA_//;s/_Metru//" all_join.txt.xls.count
selectItem.pl -h -k all_join.txt.xls.count newly_evolved_DDCP.txt |less
