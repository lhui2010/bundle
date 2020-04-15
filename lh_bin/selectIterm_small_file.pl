#!/usr/bin/perl -w


if(@ARGV <2)
{
        print "Usage: \n\t$0 list input >output\n";
        exit;
}

my ($list, $dict) = @ARGV;
my %hash;

open DICT, $dict or die;
while(<DICT>)
{
        my $t = (split /\t+/, $_)[0];
        $hash{$t} = $_;
}

open LIST, $list or $name{$list}=1;
#open IN, $ARGV[1] or die;

while(<LIST>)
{
        chomp;
       @e=split;
       print $hash{$e[0]};
}


