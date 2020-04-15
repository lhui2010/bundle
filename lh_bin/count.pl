#!/usr/bin/env perl

my $col = shift;


while(<>)
{
    chomp;
    my @e=split;
    $hash{$e[$col]}++;
}

for my $k (sort {$hash{$b} <=>$hash{$a}} keys %hash)
{
    print $k, "\t", $hash{$k}, "\n";
}
