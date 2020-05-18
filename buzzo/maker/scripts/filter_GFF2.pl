#!/usr/bin/perl -w
use strict;


my $keys = shift @ARGV;
my $GFF = shift @ARGV;

open IN, $keys or die;

my %hash;
while(<IN>)
{
    chomp;
    my @e=split;
    $hash{$e[0]} = 1;
    $hash{$e[1]} = 1;
}

close IN;

open IN, $GFF or die;

while(<IN>)
{
    my @e=split;
    my $this_id="";

    if(/mRNA\t/ or /gene\t/)
    {
        if(/ID=(.*?);/)
        {
            $this_id = $1;
        }
    }
    elsif( $e[-1]=~/Parent/)
    {
        if(/Parent=(.*?)(;|$)/)
        {
            $this_id = $1;
        }
    }
    print $_ if (exists $hash{$this_id});
}

close IN;
