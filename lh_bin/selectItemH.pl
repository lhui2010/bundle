#!/usr/bin/perl -w

# The modified version of selectItem.pl that keeps header of the to-be-selected-file
#
$source_pos = 0;
$target_pos = 0;

$keep = 0;


if(@ARGV <2)
{
        print "Usage: \n\t$0 [-k] list input >output\n";
        exit;
}

if($ARGV[0] eq "-k")
{
    $keep = 1;
    shift @ARGV;
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

$_=<>;
print;

while(<>)
{
        chomp;
        my $t = (split /\s+/, $_)[$target_pos];

        if(exists $name{$t})
        {
            if($keep)
            {
                print $name{$t},"\t",  $_ , "\n" ; # if(exists $name{$t});
            }
            else
            {
                print $_,"\n";
            }
#        delete $name{$t};
        }
}

