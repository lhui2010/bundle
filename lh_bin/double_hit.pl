#!/usr/bin/perl -w
#347     7       0       0       0       0       0       0       +       Os01t0101600-03 354     0       354     ORY_T100013263.1        673     319     673
#query 11
#target 14


open PSL1, $ARGV[0] or die;
open PSL2, $ARGV[1] or die;

while(<PSL1>)
{
	@e=split;
	$psl1{$e[9]} = $e[13];
}

while(<PSL2>)
{
	@e=split;

#if the db-qry pair have already observed in PSL1
	if(exists $psl1{$e[13]} and $psl1{$e[13]} eq $e[9])
	{
		$psl2{$e[9]} = $e[13];
	}
}

for my $key(sort keys %psl2)
{
	print $key, "\t",  $psl2{$key}, "\n";
}
