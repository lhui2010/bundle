#!/usr/bin/perl -w
use strict;

my($seq, $head, $r, $sum, $sumGC);

$/ = '>';
while(<>)
{
	chomp;
	
	my ($head, $seq) = split /\n/, $_, 2;

	next unless($head && $seq);

	$seq =~ s/\n//g;

	$r=()= $seq =~ /C|G/gi;
	
	$sumGC+=$r;
	$sum+=length($seq);

print $head, "\t", int($r*100/length($seq)), "\n";	
}

print "TOTAL\t", int($sumGC*100/$sum), "\n";
