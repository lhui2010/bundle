#!/usr/bin/perl -w

while(<>)
{
    chomp;
    my @e=split;
    if($e[2] < $e[1])
    {
        ($e[1], $e[2]) = ($e[2], $e[1]);
    }
    print (join("\t", @e), "\n");
}
