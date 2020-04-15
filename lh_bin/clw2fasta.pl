#!/usr/bin/env perl

while(<>)
{
    next if (/\*/ or /CLUST/);
    chomp;
    next if($_ eq "");
    my @e=split;
    $seq{$e[0]}.=$e[1];
}
for my $k(sort keys %seq)
{
    print ">$k\n$seq{$k}\n";
}
