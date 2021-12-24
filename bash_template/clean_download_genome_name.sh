#!/bin/bash

PREFIX=$1
longest_transcript_from_gff.py ${PREFIX}.gff ; select_ID_from_GFF.pl ${PREFIX}.gff.longest_gene ${PREFIX}.gff > ${PREFIX}.gff.longest; mv ${PREFIX}.gff ${PREFIX}.gff.orig; mv ${PREFIX}.gff.longest ${PREFIX}.gff
gff_genome_to_genes.pl ${PREFIX}.gff ${PREFIX}.fa > ${PREFIX}.cds ; cds2aa.pl ${PREFIX}.cds > ${PREFIX}.pep
add_suffix.pl ${PREFIX}.gff ${PREFIX}.cds ${PREFIX}.pep ${PREFIX}.fa
for i in *.new; do mv ${i} ${i%.new}; done

