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
	$i=0;
	for $name(sort keys %name)
	{
#		print $_ if (/$name/);
		$i=1 if (/$name/);
	}
	print $_ unless ($i);
}
