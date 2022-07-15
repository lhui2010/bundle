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
		s/$k/$k.$hash{$k}/ig;
	}
	print;
}


