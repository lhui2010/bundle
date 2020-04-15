#!/usr/bin/perl -w
open LIST, $ARGV[0] or die;
open IN, $ARGV[1] or die;

while(<LIST>)
{
	chomp;
	$name{$_} = 1;
}

while(<IN>)
{
	/ID=(.*)-..;Name/;
	$name=$1;
	@e=split;
	if(exists $name{$e[0]})
	{
		print $_;
	}
	elsif(defined $name and exists $name{$name})
	{
		print $_;
	}
}
