#!/bin/bash

ORTHO=$1
select_fasta.pl $ORTHO ~/total.pep  > ${ORTHO}.fa
mafft ${ORTHO}.fa > ${ORTHO}.fa.aln
fasttree < ${ORTHO}.fa.aln > ${ORTHO}.fa.aln.tre
