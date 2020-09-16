#!/usr/bin/perl -w

=head1 Usage:

    build_fa.pl Fasta AGP > concaternated.fa

=head1 AGP example:

    AGP_eg
    ## AGP-version 2.0
    ## AGP constructed by RaGOO
    chr01_RaGOO	1	638245	1	W	000055F|arrow_np1212	1	638245	+
    chr01_RaGOO	638246	638345	2	N	100	scaffold	yes	align_genus
    chr01_RaGOO	638346	765055	3	W	000098F|arrow_np1212	1	126710	-
    chr01_RaGOO	765056	765155	4	N	100	scaffold	yes	align_genus

=cut

die `pod2text $0` if (@ARGV == 0);

open FA, $ARGV[0] or die;
open AGP, $ARGV[1] or die;

select FA;
$/= '>';
my %hash;
while(<FA>)
{
        chomp;
        my ($head, $seq) = split /\n/, $_, 2;

        next unless($head && $seq);

        $seq =~ s/\n//g;
        $hash{$head} = $seq;
}


$/="\n";
select STDOUT;

#length of pseudo-gaps between scaffolds;
#$len_gap=100;

my %chr_seq;
my $print;
while(<AGP>)
{
##Ragoo AGP v2.0
#chr01_RaGOO	1	    638245	1	W	000055F|arrow_np1212	1	638245	+
#chr01_RaGOO	638246	638345	2	N	100	                scaffold	yes	align_genus
#chr01_RaGOO	638346	765055	3	W	000098F|arrow_np1212	1	126710	-
    next if (/^#/);

    my ($chr, $start, $end, undef, $type, $scaff_id, undef, $scaff_len, $strand)=split;

    if($type eq "W")
    {
        my $seq=$hash{$scaff_id};
        if($strand eq "-")
        {
            &rc (\$seq);
        }
        $chr_seq{$chr}.=$seq;
    }
    elsif($type eq "N")
    {
        $chr_seq{$chr}.="N"x$scaff_id;
    }
}

for my $chr( sort {substr($a, 3)<=>substr($b,3)} keys %chr_seq)
{
    print ">$chr\n";
    &formated_print($chr_seq{$chr});
}



sub formated_print
{
    $len_input=length($_[0]);
    for ($no=0; $no<$len_input; $no+=70)
    {
        print substr($_[0], $no, 70),"\n";
    }
}


sub rc
{
        my $seq_p=shift;
        if (ref($seq_p) eq 'SCALAR') { ##deal with sequence
                $$seq_p=~tr/AGCTagctMRWSYK/TCGAtcgaKYWSRM/;
                $$seq_p=reverse($$seq_p);
        }
}


