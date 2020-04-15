#!/usr/bin/perl -w

my ($len, $seq, %hash, $total, $adder, $half, $key, $fname);

$fname = $ARGV[0];
$/ = '>';
my $max = 0;
my $min = 10000;
while(<>)
{
	($len, $seq) = split /\n/, $_, 2;

	next unless (defined $seq and $seq ne "");

	$seq =~ s/\n//g;
	$seq =~ s/\s//g;

	$len = length($seq);

	$hash {$len} += $len;
	$total += $len;

	$max = $len if ($max <$len);
	$min = $len if ($min >$len);
}

print $fname, "\t", $min,"\t",$max,"\t";
$half = $total / 10;
$i=1;
$step = $half;

for $key(sort {$b <=> $a} keys %hash)
{
	$adder += $hash{$key};
	if( $adder >= $half)
	{
		$i++;
		$half+=$step;
		print  $key, "\t";
		if($i== 10){print "\n"; exit;}
	}
}
