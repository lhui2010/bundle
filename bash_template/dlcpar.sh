#!/bin/bash

sp_tree=/nfs/liuhui/bin/bundle/bash_template/dlcpar.sptree
sp_map=/nfs/liuhui/bin/bundle/bash_template/dlcpar.txt.smap
gene_tree=$1
dlcpar search -s ${sp_tree} -S ${sp_map} -D 1 -C 0.125 ${gene_tree} -I .tre  -x 1

echo "See result in ${gene_tree%.*}.dlcsearch.locus.recon"

grep dup ${gene_tree%.*}.dlcsearch.locus.recon |grep N |sort -k2,2V
