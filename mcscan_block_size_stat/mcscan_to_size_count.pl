#!/usr/bin/perl -w

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


$prefix=shift;

open SYN, "$prefix.collinearity" or die;
open BED, "$prefix.gff" or die;

#tig00000000_pilon       evm.model.tig00000000_pilon.11  218030  221012
#tig00000000_pilon       evm.model.tig00000000_pilon.110 2485412 2493675

my (%start, %end, %start_syn, %end_syn, %ctg_syn);
while(<BED>)
{
    my ($tig, $gene, $start, $end) = split;
    $start{$gene} = $start;
    $end{$gene} = $end;
}

my $gene_key_word = "pilon";
my $mark_read = 0;

my ($ref_ctg, $qry_ctg);
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
#to avoid syntenic block with gene number > 99
#eg: 322-100:
    s/-/\t/;
    my ($aln_id, $gene_id, $ref_gene, $qry_gene, $eval)=split;
    $aln_id=~s/-//;
    my $counto=0;
    my ($this_gene, $this_ctg);
    if($ref_gene =~ /$gene_key_word/)
    {
        $this_gene = $ref_gene;
        $this_ctg = $ref_ctg;
        $counto++;
    }
    if($qry_gene =~ /$gene_key_word/)
    {
        $this_gene = $qry_gene;
        $this_ctg = $qry_ctg;
        $counto++;
    }
    next if ($counto !=1);

    $start_syn{$aln_id} = $start{$this_gene} if (!exists $start_syn{$aln_id});
    $end_syn{$aln_id} = $end{$this_gene};
    $ctg_syn{$aln_id} = $this_ctg;
}

open OUT, ">$prefix.syn_block" or die;
my $sum;
for my $k(sort {$a<=>$b} keys %start_syn)
{
    $sum+=($end_syn{$k} - $start_syn{$k});
    print OUT "$k\t$ctg_syn{$k}\t$start_syn{$k}\t$end_syn{$k}\n";
}

print "$sum\n";
