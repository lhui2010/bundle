#!/usr/bin/perl
use strict;
use warnings;

my %hash;
my %position_gene;
my %position_mRNA;
my %position_CDS;
my %strand;

while(<>)
{

	my ($contig, $source, $type, $start, $end, undef, $strand, undef, $feat) = split;
	my $ID=$feat;
	$ID=~s/;.*//;
	$ID=~s/ID=//;
	($start, $end) = ($end, $start) if ($strand eq "-");

	if($type eq "gene")
	{
		next if ($ID=~/1$/);#skip pseudogene
		push @{$hash{$contig}}, $ID;
		$position_gene{$ID}="$start\t$end";
		$strand{$ID} = $strand;
	}
	elsif($type eq "mRNA")
	{
		next;
	}
	#elsif($type eq "exon" or $type eq "three_prime_UTR" or $type eq "five_prime_UTR")
	elsif($type eq "exon" )
	{
		$ID=$feat;
		$ID=~s/.*Parent=//;
		$ID=~s/;.*//;
		$ID=~s/\.\d$//;

		$position_mRNA{$ID}{"$start\t$end"} =1;
#		if($strand eq "-")
#		{
#			unshift @{$position_mRNA{$ID}}, "$start\t$end";
##			print $ID,"\n",@{$position_mRNA{$ID}}, "\n";exit;
#		}
#		else
#		{
#			push @{$position_mRNA{$ID}}, "$start\t$end";
#		}
	}
	elsif($type eq "CDS")
	{
		$ID=$feat;
		$ID=~s/.*Parent=//;
		$ID=~s/;.*//;
		$ID=~s/\.\d$//;
		$position_CDS{$ID}{"$start\t$end"} =1;
#		if($strand eq "-")
#		{
#			    unshift @{$position_CDS{$ID}}, "$start\t$end";
#		}
#		else
#		{
#			    push @{$position_CDS{$ID}}, "$start\t$end";
#		}
	}
}

sub return_start()
{
	my $input=shift;
	if($input=~/(\d+)\t/)
	{
#		print $1; exit;
		return $1;
	}
	return 0;
}

for my $k (sort keys %hash)
{
	print ">Feature $k\n";
	for my $i(0..$#{$hash{$k}})
	{
		my $gene=$hash{$k}[$i];
		print $position_gene{$gene}, "\tgene\n";
		print "\t\t\tlocus_tag\t$gene\n";
		if(exists($position_mRNA{$gene}))
		{
			my @pos1=keys %{$position_mRNA{$gene}};
			my @pos;

			if($strand{$gene} eq "-")
			{
				@pos=sort {&return_start($b) <=> &return_start($a)} @pos1;
			}
			else
			{
				@pos=sort {&return_start($a) <=> &return_start($b)} @pos1;
			}

			my $first=shift @pos;
			print $first, "\tmRNA\n";
			for my $p (@pos)
			{
				print $p, "\n";
			}
			print "\t\t\tproduct\thypothetical protein\n";
			print "\t\t\tprotein_id\tgn1|wulab|$gene.p1\n";
			print "\t\t\ttranscript_id\tgn1|wulab|$gene.t1\n";
		}
		if(exists($position_CDS{$gene}))
		{
			my @pos1=keys %{$position_CDS{$gene}};
			my @pos;

			if($strand{$gene} eq "-")
			{
				@pos=sort {&return_start($b) <=> &return_start($a)} @pos1;
			}
			else
			{
				@pos=sort {&return_start($a) <=> &return_start($b)} @pos1;
			}
			my $first=shift @pos;
			print $first, "\tCDS\n";
			for my $p (@pos)
			{
				print $p, "\n";
			}
			print "\t\t\tproduct\thypothetical protein\n";
			print "\t\t\tprotein_id\tgn1|wulab|$gene.p1\n";
			print "\t\t\ttranscript_id\tgn1|wulab|$gene.t1\n";
		}
	}
}



#C000N   maker   gene    7375    8022    .       +       .       ID=C000N0002M0;Name=C000N0002M0
#C000N   maker   mRNA    7375    8022    .       +       .       ID=C000N0002M0.1;Name=C000N0002M0;P
#C000N   maker   exon    7375    7591    .       +       .       ID=C000N0002M0.1.exon1;Parent=C000N
#
#>Features SeqID table_name
#The SeqID must match the nucleotide sequence SeqID in the corresponding .fsa file. Example Feature Table:
#
#>Feature Sc_16 Table1
#69      543    gene
#                        gene       sde3p
#69      543    CDS
#                        product SDE3P
#                        protein_id     WS1030
#
#>Feature Cont01.00055
#10      5000    gene
#                        locus_tag    GENE_NAME
#10      500     mRNA(如果是负链，按倒序，然后数字大的在前面）
#722     1555
#2548    3901
#4400    5000
#                        product     hypothetical protein
#                        protein_id  GENE_NAME-p1
#                        transcript_id   GENE_NAME-t1
#102     500     CDS
#722     1555
#2548    3901
#4400    4566
#                        product    enolase isoform A
#                        protein_id  gnl|dbname|CCC_04562A
#                        transcript_id   gnl|dbname|mrna.CCC_04562A
#
#hash{Contig}: array of genes
#hash{gene}: array of "start  end";
#						
#
#						
#protein_id应该与pep文件中的ID相同
