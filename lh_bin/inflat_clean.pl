#!/usr/bin/perl -w

while(<>)
{
  chomp;
  if(substr($_, 0, 1) eq ">")
  {
    @description=split/\s/, $_;
    $repeat_number=pop(@description);
    $desp=pop(@description);

    $seq=<>;#with \n

    for($index=1; $index<=$repeat_number; $index++)
    {
      print "$desp.$index\n$seq";
    }
  }
}
