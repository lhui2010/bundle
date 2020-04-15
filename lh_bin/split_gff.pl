#!/usr/bin/env perl
use strict;
use warnings;

my $fname=$ARGV[0];
my %hash;
while(<>)
{
    next if (/^>/ or /^#/);
    chomp;
    my @e=split;

    my $id=$e[1];

    $_.="\n";
    $hash{$id}.=$_;
}
    

for my $k (keys %hash)
{
    open OUT, ">$fname.$k.gff" or die;
    print OUT $hash{$k};
    close OUT;
}
