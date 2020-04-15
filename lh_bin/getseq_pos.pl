#!/usr/bin/perl -w

#$chr=shift @ARGV;
$fname=shift @ARGV;
$start=shift @ARGV;
$end=shift @ARGV;

open IN, $fname or die;

$seq=">";

$_=<IN>;
chomp;
$faname=$_;

while(<IN>)
{
	chomp;
	$seq.=$_;
}

print $faname,"\t",$start,"\t",$end,"\n";
&formated_print (substr($seq, $start, $end-$start+1));

sub formated_print
{
    my (@e1, $no);
    @e1=split//, $_[0];
    for($no=0; $no<@e1; $no++)
    {
        print $e1[$no];
        print "\n" if(($no+1)%70 == 0);
    }
    print "\n" unless ($no%70 == 0);
}
