#!/bin/bash

ORTHO=$1
select_fasta.pl $ORTHO ~/total.pep  > ${ORTHO}.fa
mafft ${ORTHO}.fa > ${ORTHO}.fa.aln
#fasttree < ${ORTHO}.fa.aln > ${ORTHO}.fa.aln.tre
iqtree2 -B 1000 --modelomatic -s ${ORTHO}.fa.aln
