#!/usr/bin/perl -w

my %hash;
my $pos = 0;

if(@ARGV >1)
{
    $pos = shift;
}

while(<>)
{
    my @e=split;
    if(!exists($hash{$e[$pos]}))
    {
        print;
        $hash{$e[$pos]} = 1;
    }
}
