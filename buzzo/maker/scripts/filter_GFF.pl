#!/usr/bin/perl -w
use strict;


my $keys = shift @ARGV;
my $GFF = shift @ARGV;

open IN, $keys or die;

my %hash;
while(<IN>)
{
    chomp;
    $hash{$_} = 1;
}

close IN;

open IN, $GFF or die;

while(<IN>)
{
    my @e=split;
    my $this_id="";

    if($e[-1]=~/Parent/)
    {
        if(/Parent=(.*?)(;|$)/)
        {
            $this_id = $1;
            $this_id =~s/-mRNA-\d+//;
        }
        print $_ if (exists $hash{$this_id});
    }
    else
    {
        if(/ID=(.*?);/)
        {
            $this_id = $1;
        }
        print $_ if (exists $hash{$this_id});
    }
}

close IN;
