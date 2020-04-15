#!/usr/bin/perl -w
use strict;


if(@ARGV<2)
{
    print "$0 fasta_to_be_split number_of_sub_files\n";
    exit;
}
my ($len, $seq, %hash, $total, $adder, $half, $key, $fname);

$fname = $ARGV[0];
my $pieces = $ARGV[1];

$/ = '>';
my $max = 0;
my $min = 10000;
my $total_count;

open IN, $fname or die;

while(<IN>)
{
    chomp;
    my ($id, $seq) = split /\n/, $_, 2;
    next unless (defined $seq and $seq ne "");
    $seq =~ s/\n//g;
    $seq =~ s/\s//g;

    $total += length($seq);
}

close IN;

my $dir_name="$fname._";

my $mark=1;
#if dir exists, do not rewrite
if(-d $dir_name) 
{
    $mark=0;
}
`mkdir -p $dir_name`;
my $cumulative_id;
my $cumulative_count;
my %seq;
open IN, $fname or die;
while(<IN>)
{
    chomp;
    my ($id, $seq) = split /\n/, $_, 2;

    next unless (defined $seq and $seq ne "");

    $seq =~ s/\n//g;
    $seq =~ s/\s//g;
    $id =~ s/\s+.*//;

    my $len = length($seq);

    $cumulative_count += $len;

    $cumulative_id = int((($cumulative_count-1)/$total)*$pieces)+1;

    $seq{$cumulative_id}.=">$id\n$seq\n";
}

for my $k (sort keys %seq)
{
    if($mark)
    {
        open OUT, ">$dir_name/$k.fa";
        print OUT $seq{$k};
        close OUT;
    }
    print "$dir_name/$k.fa ";
}
