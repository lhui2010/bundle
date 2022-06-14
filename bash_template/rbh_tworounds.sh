#!/bin/bash

set -euxo pipefail

#for i in  Datisca_glomerata.pep Parasponia_andersonii.pep Myrica_rubra.pep Aquilegia_coerulea.pep ; do getRBH.pl -x 24 -p mmseqs -q ${i} -t Medicago_truncatula.pep -o ${i%.pep}.Medicago_truncatula.RBH -m ${i%.pep}.Medicago_truncatula.RBH; done  

for i in *pep
do
    pushd ${i%.pep}.Medicago_truncatula.RBH/${i%.pep}.pep

    bunzip2 -qdc Medicago_truncatula.pep.rbh.bz2  |grep RBH |cut -f1 > rbh.left.id
    bunzip2 -qdc Medicago_truncatula.pep.rbh.bz2  |grep RBH |cut -f2 > rbh.right.id

    unselect_fasta.pl rbh.left.id ../../${i%.pep}.pep > ${i%.pep}.prune.pep
    unselect_fasta.pl rbh.right.id ../../Medicago_truncatula.pep > Medicago_truncatula.prune.pep

    getRBH.pl -x 24 -p mmseqs -q ${i%.pep}.prune.pep -t Medicago_truncatula.prune.pep -o ${i%.pep}.Medicago_truncatula.prune.RBH

    bunzip2 -qdc Medicago_truncatula.pep.rbh.bz2 |grep RBH |cut -f1,2 |sed "s/$/\tRBH/" > ../${i%.pep}.Medicago_truncatula.ortho.txt
    bunzip2 -qdc ${i%.pep}.Medicago_truncatula.prune.RBH/${i%.pep}.prune.pep/Medicago_truncatula.prune.pep.rbh.bz2 |grep RBH |cut -f1,2 |sed "s/$/\tsRBH/" >> ../${i%.pep}.Medicago_truncatula.ortho.txt
    popd

done
