#!/usr/bin/perl -w

print <<USAGE if (@ARGV<2);
cmp.pl by LIUHUI
Sat Apr  9 09:23:15 CST 2011

used to compare two sample (1 name per line)

usage:
cmp.pl list_A list_B >out

output:
count_sampleA_specific	count_overlap_between_A_and_B	count_sampleB_specific
USAGE

open LISTA, $ARGV[0] or die;
open LISTB, $ARGV[1] or die;

my ($count_A_sp, $count_overlap, $count_B_sp) = (0, 0, 0);

while(<LISTA>)
{
    chomp;
	$A{$_} = 1;
	$count_A_all++;
}

while(<LISTB>)
{
    chomp;
	if (exists ($A{$_}))
	{
#		$AB{$_} = 1;#overlap
		$count_overlap++;
	}
	else
	{
		#$B{$_} = 1;
		$count_B_sp++;
	}
}

$count_A_sp = $count_A_all - $count_overlap;
print $_, "\t" for ($count_A_sp, $count_overlap, $count_B_sp);
print "\n";
