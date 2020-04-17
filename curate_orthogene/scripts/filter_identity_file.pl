#!/usr/bin/perl -w

my $cutoff = 0.5;

if (@ARGV > 0)
{
    $cutoff = shift;
}

while(<>)
{
    my @e=split;

    if($e[3] ne "null" and $e[3] > $cutoff)
    {
        print;
    }
}
