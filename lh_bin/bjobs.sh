#!/bin/bash


while true
do
    date > /ds3200_1/users_root/yitingshuang/current_jobs.txt
    bjobs >> /ds3200_1/users_root/yitingshuang/current_jobs.txt
    lsload |sort -k1,1 |grep "node0[2345]\|node1[0123]" >> /ds3200_1/users_root/yitingshuang/current_jobs.txt
    sleep 1m
done

