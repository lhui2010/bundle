#!/usr/bin/perl -w

my $source_pos = 0;
my $target_pos = 0;

if(@ARGV <2)
{
        print "Usage: 
    Selece with keywords in 0-column in Query and 0-column in Ref
        $0 0 0 list input >output
    Selece with keywords in 0-column and 1-column in Query and 0-column and 1-column in Ref
        $0 0,1 0,1 list input >output
    Do not copy
        $0 0,1 0,1 list input >output\n";
        exit;
}

my $copy=1;

for my $i (0..$#ARGV)
{
    if($ARGV[$i] eq "-n")
    {
        $copy=0;
        splice @ARGV, $i, 1;

    }
}
if(int(@ARGV)>2)
{
    $source_pos = shift;
    $target_pos = shift;
}
my @source_pos =split/,/, $source_pos;
my @target_pos =split/,/, $target_pos;

open LIST, $ARGV[0] or $name{$ARGV[0]}=1;
#open IN, $ARGV[1] or die;

#if(<LIST>)
#{

while(<LIST>)
{
    chomp;
    my @e=split;
    my $t;

    for my $pos (@source_pos)
    {
        $t .=$e[$pos];
    }
    $name{$t} = $_;
#    print $t, "\n";
}

shift @ARGV;

while(<>)
{
    chomp;
    my @e = split;
    my $t = "";

    for my $pos (@target_pos)
    {
        $e[$pos] =~ s/-R.*//;
        $t .=$e[$pos];
    }
#    print $t, "\n";
#    exit;
#        my $t = (split /\s+/, $_)[$target_pos];
#    if($t eq "")
#    {
#        print STDERR $_; exit;
#    }

    if(exists $name{$t})
    {
#        print $name{$t},"\t",  $_ , "\n" ; # if(exists $name{$t});
        print $name{$t}, "\t" if $copy;
        print $_,"\n";
#        delete $name{$t};
    }
}

