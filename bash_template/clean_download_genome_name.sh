#!/bin/bash

PREFIX=$1
cp ${PREFIX}.gff ${PREFIX}.gff.orig
longest_transcript_from_gff.py ${PREFIX}.gff 
mv  ${PREFIX}.gff.longest  ${PREFIX}.gff
# gff_genome_to_genes.pl ${PREFIX}.gff ${PREFIX}.fa > ${PREFIX}.cds 
gffread -x ${PREFIX}.cds -g ${PREFIX}.fa ${PREFIX}.gff   # ${PREFIX}.cds is output
cds2aa.pl ${PREFIX}.cds > ${PREFIX}.pep
add_suffix.pl ${PREFIX}.gff ${PREFIX}.cds ${PREFIX}.pep ${PREFIX}.fa
for i in *.new; do mv ${i} ${i%.new}; done

