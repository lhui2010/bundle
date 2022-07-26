#!/bin/bash

cut -f1 $1 |sort |uniq > $1.clean
cut -f3 $1 |sed "s/.*-+-//;s/,/\n/g"|sort |uniq |grep -v "^NA$" >> $1.clean
