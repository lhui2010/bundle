#!/bin/bash

for i in $@
do
    echo $i;
    grep "^${i}" ../*.anchors |cut -f2
done
