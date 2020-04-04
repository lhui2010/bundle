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
    my @e=split;
    my $k = $e[3];
    if (exists $hash{$k})
    {
		s/$k/$hash{$k}/ig;
	}
    s/$/\t$k/;
    
	print;
}


