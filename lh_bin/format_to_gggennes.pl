#!/usr/bin/env perl
#Input
#Abrus_precatorius-NW_020874362.1_Abpre	36277011	36280618	rna-XM_027491440.1_Abpre	0	-
#Abrus_precatorius-NW_020874362.1_Abpre	36283481	36284164	rna-XM_027490886.1_Abpre	0	+
#Abrus_precatorius-NW_020874362.1_Abpre	36284218	36286767	rna-XM_027490885.1_Abpre	0	-
#Abrus_precatorius-NW_020874362.1_Abpre	36316758	36320571	rna-XM_027491343.1_Abpre	0	-
#Abrus_precatorius-NW_020874362.1_Abpre	36328499	36332607	rna-XM_027490169.1_Abpre	0	+
#Abrus_precatorius:NW_020874362.1_Abpre	36337428	36338052	rna-XM_027491829.1_Abpre	0	-
#Abrus_precatorius-NW_020874362.1_Abpre	36346870	36347566	rna-XM_027492603.1_Abpre	0	-
#Abrus_precatorius-NW_020874362.1_Abpre	36354955	36358238	rna-XM_027491939.1_Abpre	0	+
#Abrus_precatorius-NW_020874362.1_Abpre	36360194	36363497	rna-XM_027492604.1_Abpre	0	-
#Abrus_precatorius-NW_020874362.1_Abpre	36368711	36371690	rna-XM_027491779.1_Abpre	0	-
#
# Output
#  molecule gene  start    end  strand orientation
#  Genome5 <NA> 405113 407035 forward          -1
#  Genome5 genB 407035 407916 forward          -1
#  Genome5 genC 407927 408394 forward          -1
#  Genome5 genD 408387 408737 reverse          -1
#  Genome5 genE 408751 409830 forward           1
#

print join("\t", qw/molecule gene  start    end  strand orientation/), "\n";

my %strand_format;
$strand_format{"-"} = "reverse";
$strand_format{"+"} = "forward";

my %start_hash;

while(<>)
{
    next if /^--/;
    chomp;
    my $orientation = 1;
    my ($molecule, $start, $end, $gene, undef, $strand_raw) = split;
    $molecule =~ s/.*-//;
    $molecule =~ s/.*://;
    $start_hash{$molecule} ++;
    if (length($gene) > 7)
    {
        $gene = "NA";
    }
    $strand = $strand_format{$strand_raw};
    $start = $start_hash{$molecule};
    $end = $start + 0.9;
    print join("\t", ($molecule, $gene, $start, $end, $strand, $orientation)),"\n";
}
