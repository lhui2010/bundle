#!/usr/bin/perl -w

$source_pos = 0;
$target_pos = 0;

if(@ARGV <2)
{
        print "Usage: \n\t$0 list input >output\n";
        exit;
}

if(@ARGV>2)
{
    $source_pos = shift;
    $target_pos = shift;
}

open LIST, $ARGV[0] or $name{$ARGV[0]}=1;
#open IN, $ARGV[1] or die;

#if(<LIST>)
#{

while(<LIST>)
{
        chomp;
       @e=split;
       $name{$e[$source_pos]} = $_;

}

shift @ARGV;

while(<>)
{
        chomp;
        my $t = (split /\s+/, $_)[$target_pos];

        if(exists $name{$t})
        {
#        print $name{$t},"\t",  $_ , "\n" ; # if(exists $name{$t});
        print $_,"\n";
#        delete $name{$t};
        }
}

