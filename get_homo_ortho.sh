#!/bin/bash

for i in $@
do
    #echo $i;
    grep "^${i}" gene_pav.txt |less -NS
done
