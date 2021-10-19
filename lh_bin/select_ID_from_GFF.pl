#!/usr/bin/perl

#eg:
#At110g
#At012g
open LIST, $ARGV[0] or die;
#eg:
#chr1 maker gene  1  234 -  . 
#chr1 maker mRNA  1  234 -  . 
open IN, $ARGV[1] or die;
#output:
#chr1 maker mRNA 1 234 - . ID=At110g

my (%hash, @ar);
while(<LIST>)
{
	chomp;
	next if ($_ eq "");
	my $id = (split)[0];
	$hash{$id} = 1;
	push @ar, $id;
}
while(<IN>)
{
    chomp;
    my $ID = $_;
    my $Parent = $_;

    $ID=~s/.*ID=//;
    $ID=~s/;.*//;
    $Parent=~s/.*Parent=//;
    $Parent=~s/;.*//;

    if (exists $hash{$ID} or exists $hash{$Parent})
    {
        print $_,"\n";
    }
}

