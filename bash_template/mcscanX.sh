#!/bin/bash

set -euxo pipefail

# query=Senna_tora
for query in Zenia_insignis Lysidice_rhodostegia
do
    for g in `cat /nfs/liuhui/projects/phyleg/gff3_list`; 
    do 
        cd /nfs/liuhui/projects/phyleg/mcscanX
        i=${g%.gff3}

        if [ $i == $query ]
        then
            echo "do not perform self alignment this time, jumping to next"
            sleep 1s
            continue
        fi

        # Already have in first run
        # for j in cds pep gff3; do scp  -6 yitingshuang@[2400:dd07:1003:211:a94:efff:fe51:4d82]:/ds3200_1/users_root/yitingshuang/lh/projects/phyleg/mcscanX/${i}.${j} ./; done; 
        scp  -6 yitingshuang@[2400:dd07:1003:211:a94:efff:fe51:4d82]:/ds3200_1/users_root/yitingshuang/lh/projects/phyleg/mcscanX/${i}.${query}.blast ./
        bsub -J ${i}.${query} -o ${i}.${query}.log -e ${i}.${query}.err "python -m iga.evolution.ortho mcscanx ${i} ${query}"
    done
done
