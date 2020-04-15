#!/usr/bin/perl -w

#==> ../anchor_scaffold/final_result_containing_dir/scaf_chr1.1coords.filter.filter.cns.max.add_strand.agp <==
#chr01   1       84205   1       W       scaffold1024|size84205  1       84205   -
#chr01   84206   84305   2       N       100     contig  Yes

#==> Olong.representative.transcript.gff <==
#scaffold78      GLEAN   mRNA    282169  288267  0.833564        +       .       ID=Olong01m10024664.1;
#scaffold78      GLEAN   CDS     282169  282240  .       +       0       Parent=Olong01m10024664.1;

my $assemble_name = "olv1";

print "##gff-version 3\n";

my @e;
while(<>)
{
	@e=split;
	next if ($e[4] eq "N");
	$e[5]=~s/\|.*//;
	print "$e[0]\tolv1\tscaffolds\t$e[1]\t$e[2]\t.\t$e[-1]\t.\tID=$e[5];Name=$e[5];Parent=$e[0];\n"
}

print "$e[0]\tolv1\tchromosome\t1\t$e[2]\t.\t+\t.\tID=$e[0];Name=$e[0];\n"
