#!/bin/bash

DIR=/nfs/liuhui/Results_Jun21/WorkingDirectory/OrthoFinder/Results_Jul23/Orthologues/Orthologues_Medicago_truncatula

for i in $@
do
    echo $i;
    #grep "^${i}" gene_pav.txt |less -NS
    grep "^${i}" ${DIR}/gene_pav.txt |cut -f3 |sed "s/.*-+-//; s/,/\n/g" |sort |uniq |grep -v "^NA"
done
