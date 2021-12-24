#!/bin/bash
awk '$3=="CDS"' CORNE.gff |cut -f1,4,5,9 |sed 's/ID.*;Parent=//;s/;//' > CORNE.CDS.bed
sort -k1,1 -k2,2n CORNE.CDS.bed > CORNE.CDS.sort.bed
bedtools intersect -a CORNE.CDS.bed -b repeats.merge.bed -wo > CORNE.CDS.xrepeats.bed
perl -e '
my %sum;
while(<>)
{
#	/size([0-9]*)/;
	my @e=split;
	$sum{$e[0]}+=$e[1];
}

for my $k (sort keys %sum)
{
    print $k, "\t", $sum{$k}, "\n";
}'  CORNE.CDS.xrepeats.bed.len > CORNE.CDS.xrepeats.bed.len.sum
perl -ne ' my @e=split; $div = $e[1]/$e[3]; $e[4]=$div; print(join("\t", @e[0,1,2,3,4]), "\n")' CORNE.CDS.xrepeats.bed.len.sum.with_len > CORNE.CDS.xrepeats.bed.len.sum.with_len.div_result
awk '$5 >=0.5' CORNE.CDS.xrepeats.bed.len.sum.with_len.div_result > CORNE.CDS.xrepeats.bed.len.sum.with_len.div_result.overlap_larger_0.5

