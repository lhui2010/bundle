#!/usr/bin/perl -w

=head1 Usage:

    build_fa.pl Fasta AGP > concaternated.fa

=head1 AGP example:

    AGP_eg
    Pseudo-chromosome       Start   End     Scaffold        Orientation
    chr1    1       6893269 BDFN01003030    U
    chr1    6893271 9144699 BDFN01000019    U
    chr1    9144701 12478763        BDFN01002165    U
    chr1    12478765        14983719        BDFN01000179    +
    chr1    14983721        17784955        BDFN01000191    U
    chr1    17784957        18667630        BDFN01001830    U

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
$len_gap=2;

my %chr_seq;
my $print;
while(<AGP>)
{
#Pseudo-chromosome       Start   End     Scaffold        Orientation
#chr1    1       6893269 BDFN01003030    U
#chr1    6893271 9144699 BDFN01000019    U
#chr1    9144701 12478763        BDFN01002165    U
#chr1    12478765        14983719        BDFN01000179    +
#chr1    14983721        17784955        BDFN01000191    U
#chr1    17784957        18667630        BDFN01001830    U
    next if (/^Pseudo-chromosome/);

    my ($chr, $start, $end, $scaff_id, $strand)=split;

    my $seq;

    $seq=$hash{$scaff_id};

    if($strand eq "-")
    {
        &rc (\$seq);
    }
    $chr_seq{$chr}.=$seq;
    $chr_seq{$chr}.="N"x$len_gap;
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


