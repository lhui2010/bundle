#!/usr/bin/perl -w
sub transform()
{
        my $line = shift;
        my @e = split //, $line;
        return @e;
}

#[]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]<]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]f
my @total;
my @count;

my $open_arg;
if($ARGV[0] =~ /.gz$/)
{
	$open_arg = "gzip -dc $ARGV[0]|";
}
else
{
	 $open_arg = $ARGV[0];
}

open IN, "$open_arg" or die;
while(<IN>)
{
#header
#	$_;
#sequence
	my $seq = <IN>;
#+
	$_=<IN>;
	my $qual = <IN>;

#	$fq_count++;

        #next unless (defined $qual and $qual ne "");

	chomp $seq;

        my @this_qual = &transform($seq);

        for (my $i = 0; $i <@this_qual; $i ++)
        {
		if($this_qual[$i] eq "C" or $this_qual[$i] eq "G")
		{
                	$count[$i] ++;
		}
                $total[$i] ++;
        }

        #print STDERR "\r1119576 - $. -" . $./1119576 . " % OK";

}

for (my $i = 1; $i <@count; $i ++)
{
        print $i, "\t", $count[$i]/$total[$i], "\n";
}

