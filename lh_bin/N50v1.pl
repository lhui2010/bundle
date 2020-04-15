#!/usr/bin/perl -w
use strict;

my ($len, $seq, %hash, $total, $adder, $half, $key, $fname);

$fname = $ARGV[0];
$/ = '>';
my $max = 0;
my $min = 100000000;
my $total_count;

while(<>)
{
	chomp;
	my ($id, $seq) = split /\n/, $_, 2;

	next unless (defined $seq and $seq ne "");

	$seq =~ s/\n//g;
	$seq =~ s/\s//g;
	$id =~ s/\s+.*//;

	my $len = length($seq);

	$hash {$id} += $len;
	$total += $len;
	$total_count++;

	$max = $len if ($max <=$len);
	$min = $len if ($min >=$len);
}

my $num=0;
my (%Nx, %Cx);
my @indices=map { 10 * $_ } 1 .. 9;
(@indices) = (reverse @indices);

for $key(sort {$hash{$b}<=>$hash{$a}} keys %hash)
{
	$num++;
	$adder += $hash{$key};

	for my $i (@indices)
	{
		my $cutoff = $i/100;
		if( $adder/$total >= $cutoff)
		{
			unless (exists $Nx{$i})
			{
				$Nx{$i} = $hash{$key};
				$Cx{$i} = $num;
			}
		}
	}
}


#print#
print $fname,"\n";
print "Nx\tSize\tNumber\n";
for my $i (@indices)
{
		print "N$i\t", $Nx{$i}, "\t", $Cx{$i}, "\n";
}
print "Total Length\t", $total."\t".$total_count."\n";
print "Maximum Length\t", $max, "\n";
print "Minimum Length\t", $min, "\n";

