#!/bin/bash

set -euxo  pipefail

species=rice
latin="Oryza sativa"
THREADS=100


export PATH="~/lh/anaconda3/envs/repeatmodeler/bin":$PATH
##Remove after release
mkdir -p step1_repeatmask
#cat input/transcripts/*.fa >input/transcripts.fa

cd step1_repeatmask

#一定要先建立目录
#One need to add Repeat Modeler to path
mkdir -p custom_lib.out
bsub -q Q104C512G_X4  -o output.%J -e error.%J "BuildDatabase -name ref -engine ncbi ref.fa &&RepeatModeler -engine ncbi -pa $THREADS -database ref && RepeatMasker -lib ref-families.fa ref.fa -pa $THREADS -dir custom_lib.out" |awk '{print $2}' | sed 's/>//;s/<//' >>job.id

exit

mkdir -p species_lib.out
bsub -q Q104C512G_X4  -o output.%J -e error.%J "RepeatMasker ref.fa -species $species -pa $THREADS -dir species_lib.out" |awk '{print $2}' | sed 's/>//;s/<//' >>job.id

#wait  until two complete
NUM=`wc -l job.id |cut -d ' ' -f1`
while true
do
    for j in `cat job.id`
    do 
         bjobs $j >&1 |grep "not found" >>job.summary
    done

    COUNT=`wc -l job.summary |cut -d ' ' -f1`
    if [ $COUNT == $NUM ]
        then break
    else
        echo "Not completed yet..."
        sleep 1m
    fi
    rm job.summary
done

bsub <2.process_repeat.sh
