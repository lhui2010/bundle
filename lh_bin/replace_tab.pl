#!/usr/bin/perl -w
open LIST, $ARGV[0] or die;
shift @ARGV;
while(<LIST>)
{
	chomp;
	my @e=split/\t/, $_;
	$hash{$e[0]} = $e[1];
}


while(<>)
{
#    unless(/>/)
#    {
#        print;next;
#    }
	for my $k (sort keys %hash)
	{
		s/$k/$hash{$k}/ig;
	}
	chomp;
	my @e=split/\t/, $_;
    if(exists($hash{$e[1]}))
    {
        print $e[0], "\t", $hash{$e[1]};
    }
}


