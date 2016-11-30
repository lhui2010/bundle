#!/usr/bin/perl -w
use strict;

my($seq, $head, $r, $sum, $sumGC, $seq2, $LENGTH);

my (%longest, %seq, %len);

die "Usage: $0 pattern seq.fa >result.fa\nPattern eg: GVxxWELxTxxxPW\n"if(@ARGV <2);

my $pattern = shift;

#GVxxWELxTxxxPW
$pattern=~s/x/./g;

open IN, $ARGV[0] or die;
$/ = '>';
while(<IN>)
{
	chomp;
	
	my ($head, $seq) = split /\n/, $_, 2;

	next unless($head && $seq);

	$seq2 = $seq;
	$seq2 =~ s/\n//g;

    if($seq2=~/($pattern)/g)
    {
        for( my $i=1; $i <@-; $i++)
        {
            $head.="\t$-[$i]-$+[$i]";
        }
        $seq{$head} = $seq;
    }
}


for my $g (sort keys %seq)
{
    print "$g\n$seq{$g}";
}

