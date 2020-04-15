#!/usr/bin/perl -w

#count reads number from fastq file using identifier "@"

while(<>)
{
  $count_number++ if(substr($_, 0, 1) eq '@');
}

print "reads_from_fq:$count_number\n";
