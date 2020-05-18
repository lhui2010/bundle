#!/bin/bash

#for i in cn01 cn02 cn03 cn04 fat01
for i in cn01 cn02 cn03 cn04 fat01
    do
        ssh $i "rm -rf /state/partition1/tmp/${i}* && df -h |grep tmp"
        echo "$i Done.."
        done
