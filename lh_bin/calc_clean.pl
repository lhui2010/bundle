#!/usr/bin/perl -w

while(<>)
{
  chomp;
  if($_)
  {
    @lino=split /\t/, $_;
    $total_reads+=$lino[1];
    $total_bp+=$lino[0]*$lino[1];
  }
}

print "total_bp: $total_bp\n";
print "total_reads: $total_reads\n";
