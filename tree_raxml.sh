#!/bin/bash

ORTHO=$1
THREADS=$2
select_fasta.pl $ORTHO ~/total.pep  > ${ORTHO}.fa
mafft ${ORTHO}.fa > ${ORTHO}.fa.aln
#fasttree < ${ORTHO}.fa.aln > ${ORTHO}.fa.aln.tre
#iqtree2 -B 1000 --modelomatic -s ${ORTHO}.fa.aln
#iqtree2 -B 1000 -m LG+G -s ${ORTHO}.fa.aln
i=${ORTHO}.fa.aln
raxmlHPC-PTHREADS-AVX2 -f a -x  32431 -p 43213 -# 100 -m PROTGAMMAAUTO -T ${THREADS} -s ${i} -n ${i%.*}
ln -s RAxML_bipartitions.${ORTHO}.fa  ${ORTHO}.fa.aln.treefile
