#!/usr/bin/perl

open LIST, $ARGV[0] or die;
open IN, $ARGV[1] or die;
while(<LIST>)
{
	chomp;
	next if ($_ eq "");
	my $id = (split)[0];
	$hash{$id} = 1;
	push @ar, $id;
}
while(<IN>)
{
	for my $i(@ar)
	{
		if (/$i/)
		{
			print ;
			last;
		}
	}
}

