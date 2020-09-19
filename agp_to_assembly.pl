#!/usr/bin/perl -w

my $count = -1;
my $line = 0;
my @tigid;
my @list;

#my $file_len = shift;

#AGP
#Hic.fastq.gz.counts_GATC.20g10	1	132119	1	W	000083F|arrow_np1212	1	132119	-
#Hic.fastq.gz.counts_GATC.20g10	132120	132219	2	U	100	contig	yes	map
#Hic.fastq.gz.counts_GATC.20g10	132220	266408	3	W	000093F|arrow_np1212	1	134189	+
#Hic.fastq.gz.counts_GATC.20g10	266409	266508	4	U	100	contig	yes	map
#Hic.fastq.gz.counts_GATC.20g10	266509	1117534	5	W	000053F|arrow_np1212	1	851026	-
#Hic.fastq.gz.counts_GATC.20g10	1117535	1117634	6	U	100	contig	yes	map
#Hic.fastq.gz.counts_GATC.20g10	1117635	3119839	7	W	000040F|arrow_np1212	1	2002205	+
#Hic.fastq.gz.counts_GATC.20g10	3119840	3119939	8	U	100	contig	yes	map
#Hic.fastq.gz.counts_GATC.20g10	3119940	3222148	9	W	000101F|arrow_np1212	1	102209	+
#Hic.fastq.gz.counts_GATC.20g10	3222149	3222248	10	U	100	contig	yes	map
#chr01   1       4682541 1       W       000016F|arrow_np1212    1       4682541 +
#chr01   4682542 4682641 2       N       100     scaffold        yes     align_genus

my %strand;
$strand{'-'} = -1;
$strand{'+'} = 1;

my %hash;
my %len;

while(<>)
{
    next if(/^#/);
    chomp;
    my @e=split;
    next unless ($e[4] eq 'W');
    if(!exists $hash{$e[0]})
    {
        $count ++;
        $hash{$e[0]} = 1;
    }
    $line ++;
    #contig name into a list
    push @tigid, $e[5];

    #length of contig
    $len{$e[5]} = $e[-2];

    my $tmp_id = $line * $strand{$e[8]};
#    print ($line, "\n", $e[1], '\n', $strand{$e[1]}, "\n", $tmp_id, "\n");

    $list[$count] .= $tmp_id. " ";

#    >tig00002059|arrow_np1212 46 12648
#    -38 -31 28 -3 -17 25 -7 -8 12 15 -11 -2 35 -41 -36 -30 39 21 23 -34 26 29 33 -27 37 -42 -22 -24 32 20 18 19 1 5 10 16 14 9 -13 6 4
#    40

}

for my $i (0..$#tigid)
{
    print ">", $tigid[$i], "\t", $i + 1, "\t", $len{$tigid[$i]}, "\n";
}

for my $c (@list)
{
    print $c, "\n";
}
