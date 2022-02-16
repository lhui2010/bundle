#!/usr/bin/env perl

#Filter out orthologs with ten or more hits.

my (%count, %print);
while(<>)
{
    chomp;
    my @e=split;
    $count{$e[0]} += 1;
    $print{$e[0]} .= $_."\n";
}

for my $k (sort keys %count)
{
    if($count{$k} >=10)
    {
        next;
    }
    else
    {
        print $print{$k};
    }
}
