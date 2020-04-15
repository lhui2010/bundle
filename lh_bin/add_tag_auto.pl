#!/usr/bin/perl -w
$f = $ARGV[0];
$f =~s/.*\///;
$f =~ /^(...).*?_(..)/;
$f = uc($1.$2);

$tag = $f;

if (-e "$f.pep")
{
	print STDERR "please delete $f.pep before run $0";
	exit;
}
open OUT, ">$f.cds" or die;

while(<>)
{
	#if(/>/)
	if(substr($_,0,1) eq ">")
	{
		my @e=split;
		print OUT $e[0]."_$tag\n";
	}
	else
	{
		print OUT $_;
	}
}
	

