#!/usr/bin/perl -w

if(@ARGV<2)
{
	print " $0 category.txt tree.svg >tree_rename.svg\n";
	exit;
}


my $category = shift;
my $svg = shift;

my %hash;



open CAT, $category or die;
while(<CAT>)
{
	chomp; my @e=split /\t/, $_;

	$hash{$e[0]} = $e[-1];

}

open SVG, $svg or die;

while(<SVG>)
{
	for my $k(keys %hash)
	{
		s/$k/$hash{$k}/;
	}
	print;
}
