#!/usr/bin/env perl

my $qry = $ARGV[0];
my %hash;
$_=<>;
chomp;
my @f=split;
#my $qry = $f[1];
my %count;
while(<>)
{
    chomp;
    my @e=split;
    if(!exists($hash{$e[1]}))
    {
        $hash{$e[1]} = $e[0];
        $count{$e[1]} ++;
    }
    else
    {
        $hash{$e[1]} .= ",".$e[0];
        $count{$e[1]} ++;
    }
}

for my $k (sort keys %hash)
{
    print "$k\t$hash{$k}\n" if $count{$k} == 2;
}



