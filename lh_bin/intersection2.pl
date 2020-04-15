#!/usr/bin/perl -w

if (@ARGV<1)
{
	print "Usage: $0 xx yy >xx^yy\nplease input at least one file\n";
	exit;
}

$fcount = @ARGV;

#$refname = shift @ARGV;
#open REF, $refname or die;
#$current_file ++;

#while(<REF>)
#{
#	chomp;
#	@e=split;
#
#	$count{$e[0]}=$current_file;
#}

@fnames = @ARGV;
my %print;
for $fname (@fnames)
{
	open IN, $fname or die;
	
	my %tmpcount;
	
	while(<IN>)
	{
		chomp;
		next if ($_ eq "");
		@e=split;
		$tmpcount{$e[0]}=1;#$current_file if(exists $count{$e[0]});
        $print{$e[0]} .= "\t$e[1]";
		#$tmpcount{$_}=1;#$current_file if(exists $count{$e[0]});
	}

	for my $key (sort keys %tmpcount)#avoid repeated iterms in a single file
	{
		$count{$key} ++;
	}
}

for my $name_tmp(sort keys %count)
{
	if ($count{$name_tmp} >= $fcount)
    {
        print $name_tmp, $print{$name_tmp}, "\n" 
    }
}
