#!/usr/bin/env perl

while(<>)
{
    next unless (/mRNA\t/);
    $ID=$1 if (/ID=(.*?);/);
    $Name=$1 if (/Name=(.*?);/);
    $gene=$1 if (/gene=(.*?);/);
    $product=$1 if(/product=(.*?);/);

    print "$ID\t$Name\t$gene\t$product\n";
}
