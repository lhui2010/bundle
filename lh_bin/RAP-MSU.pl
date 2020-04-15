#!/usr/bin/perl -w
open INFO, "/home/SCE/liuhui/bin/RAP-MSU.txt" or die;

while(<INFO>)
{
	chomp;
	@e=split /\s+/, $_, 2;
	$info{$e[0]} = $e[1];
	
	$e[1] =~ s/\..*//;
	$info{$e[1]} = $e[0];
}

while(<>)
{
	chomp;
	@e=split;
	$e[0] =~ s/-..$//;
	$e[0] =~ s/t/g/g;
	print $e[0];
	print "\t$info{$e[0]}" if (exists $info{$e[0]});
	print "\n";
}
