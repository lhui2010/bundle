#!/usr/bin/perl -w
sub transform()
{
        my $line = shift;
        my @e = split //, $line;
        my @return = map {ord($_)-33} @e;
	warn @return;
        return @return;
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
        #$_;
#sequence
        my $seq = <IN>;
#+
        $_=<IN>;
        my $qual = <IN>;
	chomp $qual;

#        $fq_count++;
        #next unless (defined $qual and $qual ne "");

        my @this_qual = &transform($qual);

        for (my $i = 0; $i <@this_qual; $i ++)
        {
                $count[$i] ++;
                $total[$i] += $this_qual[$i];
        }

        #print STDERR "\r1119576 - $. -" . $./1119576 . " % OK";

}

for (my $i = 1; $i <@count; $i ++)
{
        print $i, "\t", $total[$i]/$count[$i], "\n";
}

