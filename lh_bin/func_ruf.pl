#!/usr/bin/perl -w
open IN,"/home/SCE/liuhui/geneloss/1_annotation_overview/annotation/3_annotation/rufipogon.function"or die;

while(<IN>)
{
	@e=split;
	$func{$e[0]} = $_;
}

while(<>)
{
	chomp;
	print $func{$_} if(exists $func{$_});
}
