#!/usr/bin/perl -w
use strict;
my %count;
my %sp;

while(<>)
{
#Group_1: b73|Zm00001d038953_P001 b73|Zm00001d011046_P001 b73|Zm00001d025188_P001 b73|Zm00001d02915
    chomp;
    my ($group_id, @rest) = split;
    $group_id=~s/://;
    $group_id=~s/.*_//;
    for my $gene (@rest)
    {   
        $gene=~s/\|.*//;
        $count{$group_id}{$gene} ++; 
        $sp{$gene}= 1;
    }   
}

#header
print "Family";
for my $sp (sort keys %sp)
{
    print "\t$sp";
}
print "\n";


for my $g (sort {$a<=>$b} keys %count)
{
    print "Group_$g";
    for my $sp (sort keys %sp)
    {   
        $count{$g}{$sp} = 0 unless (exists $count{$g}{$sp});
        print "\t$count{$g}{$sp}"; 
    }   
    print "\n";
}
