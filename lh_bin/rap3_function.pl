#!/usr/bin/perl -w
open INFO, "/home/SCE/liuhui/bin/rap3_function.data" or die;

if(@ARGV <1)
{
	print "Usage: $0 genelist >genelist_withannotation\n";
	print "This perl reads /home/SCE/liuhui/bin/rap3_function.data\n";
	print "Good Luck!\n";
	exit;
}

while(<INFO>)
{
	chomp;
	@e=split /\s+/, $_, 2;
	$info{$e[0]} = $e[1];
}
	
while(<>)
{
	chomp;
	/(Os..t.{7})/;
#tmp;
	print $_;
	print "\t",$info{$1} if (exists $info{$1});
	print "\n";
}
