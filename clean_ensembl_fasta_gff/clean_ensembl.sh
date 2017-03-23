#!/bin/bash
#A shell script used to format gene name and extract longest transcript from fasta, cds, pep, gffs of ensembl
#Dependency: BioPython, GFFparsing

#usage: edit config section -> bash clean_ensembl.sh -> success
#output example: soltu.gff3 soltu.fa soltu.pep soltu.cds soltu.gff

####config section start####
PEP=Solanum_tuberosum.SolTub_3.0.pep.all.fa
CDS=Solanum_tuberosum.SolTub_3.0.cds.all.fa
GENOME=Solanum_tuberosum.SolTub_3.0.dna_sm.toplevel.fa
GFF3=Solanum_tuberosum.SolTub_3.0.32.gff3

PREFIX=soltu

####config section end####

#TODO
#move config section to arguments


#eg. PREFIX=arath. Ath01g0120100 -> arath_Ath01g0120100. chr01 -> arath_chr01.
#prefix=$1
cur_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"



#Step1. clean gff3 format. Remove transcript: gene: prefix in ID="

cat $GFF3 |grep -v "^#" \
    |sed  "s/transcript://"  \
    |sed  "s/gene://"  \
    |sed "s/^/${PREFIX}_/" \
    |sed "s/ID=/ID=${PREFIX}_/" \
    |sed "s/Parent=/Parent=${PREFIX}_/"  \
            >$PREFIX.gff3.full

sed "s/>/>${PREFIX}_/" $GENOME >$PREFIX.fa
sed "s/>/>${PREFIX}_/" $PEP >$PREFIX.pep.full
sed "s/>/>${PREFIX}_/" $CDS >$PREFIX.cds.full


#Step2. extract longest -> PREFIX.list
python $cur_dir/longest_transcript_from_gff.py $PREFIX.gff3.full >$PREFIX.list
perl $cur_dir/extract_longest_gff.pl $PREFIX.list $PREFIX.gff3.full >$PREFIX.gff3
perl $cur_dir/select_fasta.pl $PREFIX.list $PREFIX.cds.full >$PREFIX.cds
perl $cur_dir/select_fasta.pl $PREFIX.list $PREFIX.pep.full >$PREFIX.pep

#Step3. other formats
perl $cur_dir/format_mcscan_gff_gene.pl $PREFIX.gff3 >$PREFIX.gff

