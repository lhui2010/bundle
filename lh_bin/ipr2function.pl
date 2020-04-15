#!/usr/bin/perl -w

while(<>)
{
	chomp;
	next if (/NA/ or /NULL/);
	my @e= split /\t/, $_;
	
	print $e[0], "\t";

	if(defined $e[5])
	{
		print $e[5], "\t\t\t";
	}
	
	if(defined $e[13])
	{
		print $e[13], "\n";
	}
}
