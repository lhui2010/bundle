#!/usr/bin/perl -w
$f = $ARGV[0];
$f =~s/\..*//;

$tag = $f;

while(<>)
{
	#if(/>/)
	if(substr($_,0,1) eq ">")
	{
		my @e=split;
		print $e[0]."_$tag\n";
	}
	else
	{
		print;
	}
}
	

