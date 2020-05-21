#!/usr/bin/env perl
#
if (@ARGV < 1)
{
    print 'Usage: wccount.pl  "chr*_gene_diff.tab.out.*.genblast/NS1.[34]*"\n';
    exit();
}

my $input = $ARGV[0];

my $tmp = "tmp.".rand();
system("wc -l $input > $tmp");


#my @stdout = split/\n/, $stdout;

my %sum;
open IN, $tmp or die;

while(<IN>)
{
    chomp;
    my ($count, $string)=split;
    if($string eq "total")
    {
#        print STDERR "break at ".$_."\n";
#        print STDERR "$string\n";
        next;
    }
    else
    {
        $string =~ s/.*\///;
        $sum{$string} += $count;
    }
}

for my $k (sort keys %sum)
{
    print $k, "\t", $sum{$k}, "\n";
}
system("rm $tmp");
