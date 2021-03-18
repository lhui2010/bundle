#!/usr/bin/perl -w

###gff-version|arrow_np1212	3
###sequence-region|arrow_np1212	000000F	1	12350767
#000000F|arrow_np1212	RepeatMasker	dispersed_repeat	10639	10689	230	-	.	Target=rnd-4_family-375	434	485;ID=1
#000000F|arrow_np1212	RepeatMasker	dispersed_repeat	10705	10774	478	-	.	Target=rnd-5_family-486	125	287;ID=2
#000000F|arrow_np1212	RepeatMasker	dispersed_repeat	10775	10955	578	-	.	Target=rnd-1_family-207	773	947;ID=3
#000000F|arrow_np1212	RepeatMasker	dispersed_repeat	10956	10984	478	-	.	Target=rnd-5_family-486	52	124;ID=4
#000000F|arrow_np1212	RepeatMasker	dispersed_repeat	11420	11465	583	+	.	Target=rnd-6_family-909	92	158;ID=5
#000000F|arrow_np1212	RepeatMasker	dispersed_repeat	11537	11750	583	+	.	Target=rnd-6_family-909	159	476;ID=6
#000000F|arrow_np1212	RepeatMasker	dispersed_repeat	11752	11785	466	+	.	Target=rnd-6_family-909	584	666;ID=7
#000000F|arrow_np1212	RepeatMasker	dispersed_repeat	11883	11912	466	+	.	Target=rnd-6_family-909	667	746;ID=8

while(<>)
{
    chomp;
    next if (/^#/ or $_ eq "");
    my @e=split(/\t/, $_, 9);
    my ($target, $id) = split/;/,$e[-1];
    $target =~s/Target=//;
    $target=~s/\s+/_/g;
#    print($target,"\n");
    $id =~ s/ID=//;
    my $ctg_id = $e[0];
    my $source = $e[1];
    my $repeat_type = $e[2];
    my $start = $e[3];
    my $end = $e[4];
    my $score = int($e[5]); #4.5 is not supported by maker
    my $strand = $e[6];
    my $phase = $e[7];

    ($id)=join("_", ($target, $ctg_id, $start, $end, $id));
    my $match_line = join("\t", ($ctg_id, $source, "match", $start, $end, $score,
            $strand, $phase, "ID=$id;Name=$id;Target=$target"));
    my $match_part_line = join("\t", ($ctg_id, $source, "match_part", $start, $end, $score,
            $strand, $phase, "ID=$id.part;Parent=$id"));
    print($match_line, "\n", $match_part_line, "\n");
}

