#!/usr/bin/perl

$_=<>;
my @e=split;
for (my $i=1; $i <=int(@e); $i++)
{
    print $i, "\t";
}
print "\n";
print;

my $count=5;
while(<>)
{
    if($count--<0)
    {
        last;
    }
    else
    {
        print;
    }
}

