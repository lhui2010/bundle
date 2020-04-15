#!/usr/bin/perl -w
#usage: find_unmapped.pl xx.bowtie xx.fa
#find unmapped reads by comparing bowtie result and clean.fa
#then print it 
#pipe: read reads' names from .bowtie to array, sort, then compare
#15:46 2010/8/4

if(!@ARGV)
{
print "find unmapped reads by comparing bowtie result and clean.fa\nusage: find_unmapped.pl xx.bowtie xx.fa\n";
exit;
}

$file_bowtie=$ARGV[0];
$file_fa=$ARGV[1];

open BOWTIE, $file_bowtie||die"can't open $file_bowtie\n";
open FA, $file_fa||die"can't open $file_fa\n";



while(<BOWTIE>)
{
@e=split/\s+/, $_;
push(@name, $e[0]);
}

@sorted_name=sort @name;


$line_name=<FA>;
$name_F=substr($line_name, 1, 8);


for $name_B(@sorted_name)
{
while ($name_F lt $name_B && $name_F ne "")
{
$line_seq=<FA>;
print "$line_name$line_seq";

$line_name=<FA>;
$name_F=substr($line_name, 1, 8);
}
$line_name=<FA>;#seq's name
$line_name=<FA>;

exit unless($line_name);

$name_F=substr($line_name, 1, 8);
}

