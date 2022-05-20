#!/usr/bin/perl -w

# Update 01.04
# Accept input like selectItem.pl 0 0,1 table target > table_in_target

my $source_pos = 0;
my $target_pos = 0;
my $keepheader = 0;
my $keep = 0;


if(@ARGV <2)
{
        print "Usage: \n\t$0 [-h][-k] list input >output\n";
        exit;
}

if($ARGV[0] eq "-h")
{
    $keepheader = 1;
    shift @ARGV;
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

my @target_pos = split/,/, $target_pos;

if($keepheader)
{
    $_=<>;
    print;
}

while(<>)
{
        chomp;
        my @t;
        for my $t (@target_pos)
        {
            #push @t, (split /\s+/, $_)[$target_pos];
            push @t, (split /\s+/, $_)[$t];
        }

        my $flag = 1;

        for my $t(@t)
        {
            if(!exists $name{$t})
            {
                $flag = 0;
            }
        }

        #if(exists $name{$t})
        if($flag)
        {
            if($keep)
            {
                print $name{$t[0]},"\t",  $_ , "\n" ; # if(exists $name{$t});
            }
            else
            {
                print $_,"\n";
            }
#        delete $name{$t};
        }
}

