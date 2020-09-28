#!/bin/bash
#BSUB -J process_repeat     # job name
#BSUB -n 1                   # number of tasks in job
#BSUB -q Q104C512G_X4              # queue
#BSUB -e errors.%J     # error file name in which %J is replaced by the job ID
#BSUB -o output.%J     # output file name in which %J is replaced by the job ID

set -euxo pipefail 

species=cucumber

cd step1_repeatmask
mkdir -p Full_mask

gunzip *lib.out/*.cat.gz
cat *lib.out/*.cat >Full_mask/full_mask.cat
cd Full_mask
#ProcessRepeats -species "Cucumis sativus" full_mask.cat
ProcessRepeats -species ${species} full_mask.cat

#create GFF3
#remove after release, use nohup to run, can't use sge
rmOutToGFF3.pl full_mask.out > full_mask.gff3

# isolate complex repeats
#remove after release, use nohup to run, can't use sge
grep -v -e "Satellite" -e ")n" -e "-rich" full_mask.gff3 \
    > full_mask.complex.gff3
# reformat to work with MAKER
  cat full_mask.complex.gff3 | \
    perl -ane '$id; if(!/^\#/){@F = split(/\t/, $_); chomp $F[-1];$id++; $F[-1] .= "\;ID=$id"; $_ = join("\t", @F)."\n"} print $_' \
            > full_mask.complex.reformat.gff3

echo "Repeat GFF file is located in "
echo "$PWD/full_mask.complex.reformat.gff3"

