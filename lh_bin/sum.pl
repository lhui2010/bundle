#!/usr/bin/perl -w

my $column = shift;

#$_=<>;
while(<>)
{
#	/size([0-9]*)/;
	my @e=split;
	$sum+=$e[$column];
}
print $sum, "\n";
