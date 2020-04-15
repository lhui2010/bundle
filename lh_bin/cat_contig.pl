#!/home/SCE/liuhui/bin/perl/bin/perl -w

if(@ARGV <2)
{
	warn "\n\nError: Insufficient parameter
Program:\n\t$0
Usage:\n\t$0 *.contig >merged.contig
\n\n";
	exit;
}

@fnames = @ARGV;

for my $fname (@fnames)
{
	open IN, $fname or die;

	my ($spname, undef) = split /\./, $fname, 2;

	while(<IN>)
	{
		if(/>/)
		{
			print ">$spname"."-";
			print substr $_, 1;
		}
		else
		{
			print $_;
		}
	}
}
