#!/usr/bin/perl -w

#Function: Extract orthologous genes identified using MCScanX
#Author: lhui2010@gmail.com
################ Parameters ###############
## MATCH_SCORE: 50
## MATCH_SIZE: 5
## GAP_PENALTY: -1
## OVERLAP_WINDOW: 5
## E_VALUE: 1e-05
## MAX GAPS: 25
################ Statistics ###############
## Number of collinear genes: 26206, Percentage: 50.47
## Number of all genes: 51922
###########################################
### Alignment 0: score=389.0 e_value=3.3e-20 N=10 1&1 plus
#  0-  0:	Solyc01g008380.1.1	Solyc01g096200.2.1	  8e-19
#  0-  1:	Solyc01g008480.2.1	Solyc01g096210.2.1	  7e-17
#  0-  2:	Solyc01g008640.2.1	Solyc01g096280.1.1	  3e-60
#  0-  3:	Solyc01g008740.1.1	Solyc01g096350.2.1	 2e-103
#  0-  4:	Solyc01g008850.2.1	Solyc01g096490.2.1	  4e-10
#  0-  5:	Solyc01g008860.1.1	Solyc01g096690.1.1	  9e-46
#  0-  6:	Solyc01g008960.2.1	Solyc01g096750.1.1	      0
#  0-  7:	Solyc01g009170.2.1	Solyc01g096810.2.1	      0

my $gene_key_word = "pilon";

my $prefix=shift;

open SYN, "$prefix.collinearity" or die;

my $mark_read = 0;

open OUT1, ">$prefix.$gene_key_word" or die;
open OUT2, ">$prefix.cross" or die;
open OUT3, ">$prefix.ref" or die;


while(<SYN>)
{
    if (/Alignment/)
    {
        $mark_read =1 ;
        my @e=split;
        $ref_ctg = $e[6];
        $ref_ctg =~s/\&.*//;
        $qry_ctg = $e[6];
        $qry_ctg =~s/.*\&//;
        next;
    }
    next unless ($mark_read);
    next if (/Alignment/);

#to avoid syntenic block with gene number > 99
#eg: 322-100:
    s/-/\t/;

#only genes extract with gene key word
    my ($aln_id, $gene_id, $ref_gene, $qry_gene, $eval)=split;
    $aln_id=~s/-//;
    my $counto=0;
    my ($this_gene, $this_ctg);
    if($ref_gene =~ /$gene_key_word/ and $qry_gene =~ /$gene_key_word/)
    {
        print OUT1 $ref_gene,"\t",$qry_gene,"\n";
    }
    elsif($ref_gene =~ /$gene_key_word/ or $qry_gene =~ /$gene_key_word/)
    {
        print OUT2 $ref_gene,"\t",$qry_gene,"\n";
    }
    else
    {
        print  OUT3 $ref_gene,"\t",$qry_gene,"\n";
    }
}
