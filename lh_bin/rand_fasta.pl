#!/usr/bin/perl -w

#打印出随机的序列，默认是选百分之一

$max_seq_for_print = 100;

$/ = '>';
<>;
while(<>)
{
	($len, $seq) = split /\n/, $_, 2;
	$buffer{$len} = $seq;
}


my $count = 0;
for my $print (keys %buffer)
{
	print ">",$print,"\n", $buffer{$print};
	$count ++;
	last if( $count > $max_seq_for_print);
}
