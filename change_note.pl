#!/usr/bin/perl

if(@ARGV <2)
{
    print "$0 -0.5 1 3 5 5 4 6 1\n";
    exit;
}

my ($change, @notes) = @ARGV;

my %decode=(
    1=>1,
    1.5=>2,
    2=>3,
    2.5=>4,
    3=>5,
    4=>6,
    4.5=>7,
    5=>8,
    5.5=>9,
    6=>10,
    6.5=>11,
    7=>12,
    );

my %encode = reverse %decode;

$change *=2;

for my $note(@notes)
{
    my $new = $decode{$note} + $change;
    $new +=12 if ($new <=0);
    print $encode{$new} , " ";
}
