#!/usr/bin/perl -w



if(@ARGV <2)
{
        print "Usage: \n\t$0 list input >output\n";
        exit;
}


#Order of column of list in input
my $index = 0;




open LIST, $ARGV[0];
#open IN, $ARGV[1] or die;

while(<LIST>)
{
      chomp;
      next if ($_ eq "");

      my @keyword = split;
      $transcript{$keyword[0]} = 1;
      $gene{$keyword[1]} = 1;
}





shift @ARGV;

#print header
#$_=<>;
#print;

while(<>)
{
    if(/\tchromosome/)
    {
        print;
    }
    elsif(/\tgene/)
    {
        /ID=(.*?);/;
        print if (exists $gene{$1});
    }
    elsif(/\ttranscript/)
    {
        /ID=(.*?);/;
        print if (exists $transcript{$1});
    }
    else
    {
        if(/Parent=(.*?);/)
        {
            print if(exists $transcript{$1});
        }
    }
}
