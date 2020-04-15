#!/usr/bin/perl -w
#use strict;
use warnings;

warn "$0 a/b coverage(two-way) .psl >result\n" if (@ARGV<2);

my $coverage = 0.3;
my $expression = "a";

$expression = shift @ARGV; 
$coverage = shift @ARGV; 

my $fname = $ARGV[0];

open OUT , ">".$fname.".$expression$coverage" or die;

$coverage /=100 if ($coverage >1);

my($match, $qname, $qsize, %below, %bsub, %asub, %above, @e);

while(<>)
{
	chomp;
	if($_ eq "" or /psLayout/ or /match/ or /-----/)
	{
		next;
	}
	@e=split;

	next if( $e[9] eq $e[13]);
	
	$match = $e[0];
	$qname = $e[9];
	$qsize = $e[10];
	$tname = $e[13];
	$tsize = $e[14];

	if ($match/$qsize < $coverage and $match/$tsize < $coverage )#and $match/$tsize < $coverage)
	{
		$below{$qname} = $match/$qsize ;
		$bsub{$qname} .= $tname;
		$bseq{$qname} = $_;
	}
	if($match/$qsize >= $coverage or $match/$tsize >= $coverage)# and $match/$tsize >= $coverage)
	#print STDERR "$match\t$tsize\t$coverage\n";
	#if($match/$tsize >= $coverage)
	{
		$above{$qname} = $match/$qsize ;
		$asub{$qname} .= "$tname\t";
		$aseq{$qname} = $_;

		#print STDERR "$qname\t$aseq{$qname}\n";
	}
}

if($expression eq "b")
{
	#print $_, "\t", $bsub{$_}, "\n" for (sort keys %below);
	for $key(sort keys %bseq)
	{
		print OUT $bseq{$key}, "\n";
	}
}
elsif($expression eq "a")
{
#	print $_, "\t", $asub{$_}, "\n" for (sort keys %above);
	for $key(sort keys %aseq)
	{
		print OUT $aseq{$key}, "\n";
	}
}


