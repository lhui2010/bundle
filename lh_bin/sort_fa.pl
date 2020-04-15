#!/usr/bin/env perl

$/='>';
$_=<>;
while (<>) 
{
	chomp;
	($h,$s)=split /\n/, $_, 2; 
	$seqs{$h}=$s;
}
foreach $header (sort keys %seqs) 
{
	print ">".$header."\n".$seqs{$header}
}
