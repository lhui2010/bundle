#!/bin/bash

sp_tree=$BD/bash_template/dlcpar.sptree
sp_map=$BD/bash_template/dlcpar.txt.smap
gene_tree=$1
#dlcpar search -s ${sp_tree} -S ${sp_map} -D 1 -C 0.125 ${gene_tree} -I .tre  -x 1
dlcpar dp -s ${sp_tree} -S ${sp_map} ${gene_tree} --output_format 3t

echo "See result in ${gene_tree%.*}.dlcdp.locus.recon"

#grep dup ${gene_tree%.*}.dlcsearch.locus.recon |grep N |sort -k2,2V
