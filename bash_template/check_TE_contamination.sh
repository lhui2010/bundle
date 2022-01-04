#!/bin/bash

# Input: 
#     1. cds
#     2. repeat_mask.out
#     3. gff

# If report perl lib not found 
# Do the following:
# cd $CONDA_PREFIX/lib/site_perl/5.26.2/
# for i in $CONDA_PREFIX/share/RepeatMasker/*.pm; do ln -s $i; done
rmOutToGFF3.pl Senna_tora.fa.out > Senna_tora.fa.gff 
awk '$3=="CDS"' Senna_tora.gff | gff2bed.pl > Senna_tora.bed
bedtools coverage -a Senna_tora.bed -b Senna_tora.fa.gff > Senna_tora.bed.xrepeat.coverage
cut -f4,8 Senna_tora.bed.xrepeat.coverage |sed "s/CDS.*\.//" > Senna_tora.bed.xrepeat.coverage.extract
sumdict.pl Senna_tora.bed.xrepeat.coverage.extract  > Senna_tora.bed.xrepeat.coverage.extract.sum
samtools faidx Senna_tora.cds
cut -f1,2 Senna_tora.cds.fai > Senna_tora.cds.fai.size
perl frac.pl Senna_tora.bed.xrepeat.coverage.extract.sum Senna_tora.cds.fai.size > Senna_tora.cds.fai.size.repeat_frac
awk '$3<0.5' Senna_tora.cds.fai.size.repeat_frac |cut -f1 > Senna_tora.cds.clean.txt
