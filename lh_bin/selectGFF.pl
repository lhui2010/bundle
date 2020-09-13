#!/usr/bin/perl -w
open LIST, $ARGV[0] or die;
open IN, $ARGV[1] or die;

while(<LIST>)
{
	chomp;
	$name{$_} = 1;
}

while(<IN>)
{
    my $name = "";
	if(/ID=(.*);Name/ or /Parent=(.*)/)
    {
        $name=$1;
    }
	@e=split;
	if(exists $name{$e[0]})
	{
		print $_;
	}
	elsif(defined $name and exists $name{$name})
	{
		print $_;
	}
}
