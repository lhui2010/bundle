#!/bin/bash

for i in $@
do
    echo $i;
    #grep "^${i}" ../*.anchors |cut -f2
    grep "^${i}" /nfs/liuhui/projects/buzzo/quota_align/*.anchors |cut -f2
done
