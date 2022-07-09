#!/usr/bin/perl
while(<>)
{
    if(/#/)
    {
        print;
    }
    else
    {
        chomp;
        my @e=split; 
        ($e[0], $e[1]) = ($e[1], $e[0]);
        $_=join("\t", @e);
        $_.="\n";
        print;
    }
}
