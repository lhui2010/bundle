#!/usr/bin/perl -w

$depth_threshold = 5;

while(<>)
{
    next if (/^#/);
    chomp;
    my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, $sample1) = split;
    if($ALT=~/,/)
    {
#GT:AD:DP:GQ:PL
###FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
        my ($GT,$depth, $DP, $GQ, $PL) = split /:/, $sample1;
        my @alt_depth = split /,/, $depth;
        my $is_ok =1;
        shift @alt_depth if ($alt_depth[0] ==0);
        for my $i (@alt_depth)
        {
            if($i < $depth_threshold)
            {
                $is_ok=0;
            }
        }
        print $_, "\n" if ($is_ok);
    }
    else
    {
        next;
    }
}
