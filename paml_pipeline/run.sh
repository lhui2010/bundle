#!/bin/bash

#lhui2010@gmail.com

#############################
# Things you need provide:  #
#                           #
# Argument:                 #
#    $group_name            #
# Files:                    #
#    PEP/$group_name.pep    #
#    CDS/$group_name.cds    #
#
# Softwares:                #
#    T_Coffee               #
#    PAL2NAL                #
#    FastTree               #
#    clustalw2fasta(meme)   #
#    codeml (PAML)          #
#############################

g=$1

#echo $g
#exit
cd $PWD

mkdir work/$g

#input pep and cds
cd work/$g/
ln -s ../../PEP/${g}.pep
ln -s ../../CDS/${g}.cds

t_coffee ${g}.pep -method mafftgins_msa,muscle_msa,kalign_msa,t_coffee_msa > ${g}.pep.aln
pal2nal.pl ${g}.pep.aln ${g}.cds   >${g}.paml_aln
clustalw2fasta ${g}.paml_aln >${g}.paml_aln.fa
fa2phy.pl ${g}.paml_aln.fa >${g}.paml_aln.phy #OK
FastTreeMP <${g}.paml_aln.fa >${g}.paml_tree #wrong format, use phylip
#remove bootstrap value
sed -ie 's/)[0-9].[0-9]\+:/):/g' ${g}.paml_tree #remove bootstrap value
#add #1 to genes of interests
sed -ie 's/\(Cau|C[0-9]\+[A-Z][0-9]\+[A-Z][0-9].[0-9]\)/&#1/'g ${g}.paml_tree

echo "seqfile = ${g}.paml_aln.phy
treefile = ${g}.paml_tree
outfile = alte.mlc

noisy = 3
verbose = 1
runmode = 0

seqtype = 1
CodonFreq = 2
clock = 0
model = 2

NSsites = 2
icode = 0

fix_kappa = 0
kappa = 2.5
fix_omega = 0
omega = 1.5

fix_alpha = 1
alpha = .0
Malpha = 0
ncatG = 4

getSE = 0
RateAncestor =0

fix_blength = -1
method = 0
">alte.ctl

echo "seqfile = ${g}.paml_aln.phy
treefile = ${g}.paml_tree
outfile = null.mlc

noisy = 3
verbose = 1
runmode = 0

seqtype = 1
CodonFreq = 2
clock = 0
model = 2

NSsites = 2
icode = 0

fix_kappa = 0
kappa = 2.5
fix_omega = 1
omega = 1

fix_alpha = 1
alpha = .0
Malpha = 0
ncatG = 4

getSE = 0
RateAncestor =0

fix_blength = -1
method = 0
">null.ctl

codeml null.ctl
codeml alte.ctl

