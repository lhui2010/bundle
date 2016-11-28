#!/usr/bin/perl -w

#my $column = shift;
my $end = 3;
my $start = 2;

while(<>)
{
#	/size([0-9]*)/;
	my @e=split;
#	$sum+=$e[$column];
	$sum+=$e[$end] - $e[$start];
}
print $sum, "\n";
