#!/usr/bin/perl -w 
while(<>)
{
	chomp;
	my @e=split /\t/, $_;

	print $_, "\n" if ($e[4] >0.5)# and $e[5] < 0.1);
}
