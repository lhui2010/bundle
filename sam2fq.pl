#!/usr/bin/perl

while(<>)
{
	chomp;
	my @e=split;
	my $id=$e[0];
	my $nuc = $e[9];
	my $qual = $e[10];

	print "\@$id\n$nuc\n+$id\n$qual\n";
}
