#!/bin/bash
#BSUB -J kmer      # job name
#BSUB -n 1                   # number of tasks in job
#BSUB -q Q104C512G_X4              # queue
#BSUB -e errors.%J     # error file name in which %J is replaced by the job ID
#BSUB -o output.%J     # output file name in which %J is replaced by the job ID

set -euo pipefail 

ROOT=$PWD
WORKDIR=${ROOT}/workdir0805
SCRIPTDIR=${ROOT}/bin

mkdir -p $WORKDIR && cd ${WORKDIR}

#export LD_LIBRARY_PATH=""

echo -e "Species\tGenomeSize\tHeterozygosity\tRepeat%"

for i in `cat ${ROOT}/sample_list`
do
#    ls ${ROOT}/input/${i}*/*gz > ${i}.fq.lst
#    ${SCRIPTDIR}/gce-1.0.0/kmerfreq/kmer_freq_hash/kmer_freq_hash -t 64 -k 23 -l ${i}.fq.lst -p ${i} 2>${i}.kmerfreq.log
    UNIQKMERNUM=`tail -n 11 ${i}.kmerfreq.log  |head -1 |awk '{print $2}'`
    DEPTH=`tail -n 11 ${i}.kmerfreq.log  |head -1 |awk '{print $5}'`
    #${SCRIPTDIR}/gce-1.0.0/gce -f ${i}.freq.stat -g ${UNIQKMERNUM} -H 1 -c ${DEPTH} -b 1 >${i}.gce.out 2>${i}.gce.err
# Genome Size
#    echo -en "${i}\tGenome Size:\t"; tail -n 2 ${i}.gce.err |head -1 |awk '{print $6}'
# Heterozygosity 
#    echo -en "\tHeterozygosity:\t"; tail -n 2 ${i}.gce.err |head -1 |awk '{print $7/(2-$7)/23}' 
# Repeat Rate
#    echo -en "\tRepeat%:\t"; tail -n 2 ${i}.gce.err |head -1 |awk '{print 1-$9-$10}' 
    GenomeSize=`tail -n 2 ${i}.gce.err |head -1 |awk '{print $6}'`
    Heterozygosity=`tail -n 2 ${i}.gce.err |head -1 |awk '{print $7/(2-$7)/23}'`
    Repeat=`tail -n 2 ${i}.gce.err |head -1 |awk '{print 1-$9-$10}'`
    echo -e "${i}\t${GenomeSize}\t${Heterozygosity}\t${Repeat}"
done

exit
#./bin/gce-1.0.0/kmerfreq/kmer_freq_hash/kmer_freq_hash -t 64 -k 23 -l fq.lst 2>kmerfreq.log
#Manual check kmerfreq.log and get total kmer number, here is 45969582290
#./bin/gce-1.0.0/gce -f output.freq.stat  -g 45969582290 -H 1 -c 126 -b 1 >gce.out 2>gce.err
#./bin/gce-1.0.0/gce -f 23mer.freq  -H 1 -c 126

#Calculate Heterozygosity
#-H means calculate heterozygosity rate H, which is H = a[1/2] / (2-a[1/2]) / 23, where 23 is kmer length

