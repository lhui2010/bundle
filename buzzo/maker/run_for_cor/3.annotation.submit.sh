#!/bin/bash

set -euxo pipefail

DIR=$PWD

#cp -r /ds3200_1/users_root/yitingshuang/lh/bin/maker3/exe/augustus/config /tmp/lh_config
export AUGUSTUS_CONFIG_PATH=/tmp/lh_config

#remove after release
#cat input/peps/*.fa >input/peps.fa
#cat input/transcripts/*.fa >input/transcripts.fa

#mv step2_annotation step2_annotation.rm
mkdir -p step2_annotation
##
cd step2_annotation
#

##generating CTL files

OPT=round1_maker_opts.ctl

for i in `cat fail_node10.txt.id`
    do
        cd $DIR/step2_annotation/round1
        rm -rf maker.$i
        mkdir -p maker.$i
        cd maker.$i
        ln -s ../../fasta_partition/$i.fa
        cp ../../${OPT} ./
        cp ${DIR}/maker_exe.ctl ./
        cp ${DIR}/maker_bopts.ctl ./
        echo "genome=$PWD/$i.fa #genome sequence (fasta file or fasta embeded in GFF3 file)" >>$OPT
#3.01.02
        #bsub -q Q104C512G_X4 -J R1mk$i -o maker.$i.log -e maker.$i.err " export AUGUSTUS_CONFIG_PATH=/tmp/lh_config && /ds3200_1/users_root/yitingshuang/lh/bin/maker/bin/maker $OPT maker_bopts.ctl maker_exe.ctl"
        bsub -q Q104C512G_X4  -J mk$i -o output.%J -e error.%J "export AUGUSTUS_CONFIG_PATH=/tmp/lh_config && /ds3200_1/users_root/yitingshuang/lh/bin/maker_nompi/bin/maker *ctl >maker.${i}.log 2> maker.${i}.err "
        sleep 3s
#        echo "/lustre/local/packages/Maker/maker-3.01.02/bin/maker $OPT maker_bopts.ctl maker_exe.ctl" >run.sh
    done
