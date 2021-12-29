#  ==> Medicago_truncatula.Medicago_truncatula.tdarray <==
#  Anchor1	Tandem_array1	Anchor2	Tandem_array2
#  mRNA:MtrunA17CPg0493221_Metru	mRNA:MtrunA17CPg0493221_Metru,mRNA:MtrunA17CPg0493231_Metru	mRNA:MtrunA17Chr4g0015351_Metru	mRNA:MtrunA17Chr4g0015351_Metru
#  mRNA:MtrunA17CPg0493221_Metru	mRNA:MtrunA17CPg0493221_Metru,mRNA:MtrunA17CPg0493231_Metru	mRNA:MtrunA17Chr4g0016511_Metru	mRNA:MtrunA17Chr4g0016511_Metru,mRNA:MtrunA17Chr4g0016521_Metru
#  mRNA:MtrunA17CPg0493221_Metru	mRNA:MtrunA17CPg0493221_Metru,mRNA:MtrunA17CPg0493231_Metru	mRNA:MtrunA17Chr7g0229251_Metru	mRNA:MtrunA17Chr7g0229251_Metru
#  mRNA:MtrunA17Chr0c01g0489031_Metru	mRNA:MtrunA17Chr0c01g0489031_Metru	mRNA:MtrunA17Chr2g0306241_Metru	mRNA:MtrunA17Chr2g0306241_Metru,mRNA:MtrunA17Chr2g0306251_Metru

# ==> Senna_tora.Medicago_truncatula.collinearity <==
#  ############### Parameters ###############
#  # MATCH_SCORE: 50
#  # ...
#  ## Alignment 0: score=384.0 e_value=2e-14 N=8 3___fragment_2___debris__unscaffolded_Setor&MtrunA17Chr7_Metru plus
#    0-  0:	Sto00g444800_Setor	mRNA:MtrunA17Chr7g0276201_Metru	  1e-09
#    0-  1:	Sto00g444830_Setor	mRNA:MtrunA17Chr7g0276281_Metru	  4e-54
  

die "$0 tdarray collinearity " if(@ARGV < 2);

my $tdarray_file = shift;
my $collinearity_file = shift;

`col2anc.pl $collinearity_file > $collinearity_file.anchors`;
`grep -v "#" $collinearity_file.anchors |awk '{print \$1"\t"\$2}' > $collinearity_file.anchors.ortho`;

my %remove_dict;
open $IN, $tdarray_file;
while(<IN>)
{
    chomp;
    my ($key_ref, $bak_ref, $key_qry, $bak_qry)=split;
    my (@bak_ref, @bak_qry);
    if($bak_ref =~/,/)
    {
        @bak_ref = split/,/, $bak_ref;
    }
    else
    {
        @bak_ref = ($bak_ref);
    }

    if($bak_qry =~/,/)
    {
        @bak_qry = split/,/, $bak_qry;
    }
    else
    {
        @bak_qry = ($bak_qry);
    }

    for my $bq (@bak_qry)
    {
        $remove_dict{$key_ref}{$bq} = 1;
        $remove_dict{$bq}{$key_ref} = 1;
        for my $br (@bak_ref)
        {
            $remove_dict{$br}{$bq} = 1;
            $remove_dict{$bq}{$br} = 1;
        }
    }
    for my $br (@bak_ref)
    {
        $remove_dict{$key_qry}{$br} = 1;
        $remove_dict{$br}{$key_qry} = 1;
    }
}
close IN;

open IN, "$collinearity_file.anchors.ortho" or die;
while(<IN>)
{
    chomp;
    my @e=split;
    if(exists $remove_dict{$e[0]}{$e[1]})
    {
        next;
    }
    print $e[0], "-", $e[1], "\n";
}


