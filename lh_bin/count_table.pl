#!/usr/bin/perl
while(<>)
{
	chomp;
	my @e=split;
	my ($start, $end) = ($e[0], $e[1]);
	my $minus = abs($end-$start);
	$sum+=$minus;
}
print $sum, "\n";
