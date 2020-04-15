#!/usr/bin/perl -w

$/='>';
while(<>)
{
	chomp;
	next if ($_ eq "");

	my ($desc, $fa) =split (/\n/, $_, 2);

	$fa=~s/\n//g;
	$fa=~s/\r//g;
	$desc=~s/\r//g;

	print ">$desc\n";
	&rc(\$fa);
	print $fa,"\n";
}

sub rc
{
	my $seq_p=shift;
	if (ref($seq_p) eq 'SCALAR') { ##deal with sequence
	    $$seq_p=~tr/AGCTagctMRWSYK/TCGAtcgaKYWSRM/;
	    $$seq_p=reverse($$seq_p);
	}
}
