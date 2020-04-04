#!/usr/bin/perl -w

$source_pos = 0;
$target_pos = 0;

if(@ARGV<2)
{
	        print "$0 special total >total_without_special\n";
		        exit;
}
elsif(@ARGV== 4)
{
    $source_pos = shift;
    $target_pos = shift;
}

open LIST, $ARGV[0] or die;
open IN, $ARGV[1] or die;

while(<LIST>)
{
	        chomp;
	        next if ($_ eq "");
	        $name{(split)[$source_pos]} = 1;
}

#$target_pos = 0;
while(<IN>)
{
	s/>//;
    my @e=split;
        unless(exists $name{$e[$target_pos]})
        {
                print $_ ;#if ($mark);
#                $name{$e[$target_pos]} = 1
        }
}

