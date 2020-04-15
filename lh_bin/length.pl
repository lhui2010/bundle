#!/usr/bin/perl -w
use strict;

my ($len, $seq, %hash, $name, $adder, $half, $key, $fname);

$fname = $ARGV[0];
$/ = '>';
while(<>)
{
        my $tmp;
        ($tmp, $seq) = split /\n/, $_, 2;

        next unless (defined $seq and $seq ne "");

        $name = (split /\s+/, $tmp) [0];

        $seq =~ s/\n//g;
        $seq =~ s/\s//g;

        $hash{$name} = length($seq);
}


for $key(sort keys %hash)
{
        print $key, "\t", $hash{$key}, "\n";
}
