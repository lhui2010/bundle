#!/usr/bin/env perl
# add  suffix to Arabidopsis.gff and Arabidopsis.fa
# Input: >chr1
# Output: >chr1_artha
# Input: chr1  maker ... ID=AT1G0001
# Output: chr1_artha maker ... ID=AT1G0001_artha
# Suffix is defined in line 17
# Usage: perl add_suffix.pl Arabidopsis.gff Arabidopsis.fa Arabidopsis.cds Arabidopsis.pep
# Input file are formated with species name and format extension.
# New file will be written to Arabidopsis.gff.new Arabidopsis.fa.new, etc

if (@ARGV < 1)
{
    die "add_suffix.pl Vitis_vinifera.*"
}

@files=@ARGV;


#TODO: 把特殊符号给替换了，比如 | :，不做的话，一些软件比如KaKs和唐海宝的流程都会改成下划线，造成基因名称无法对称的情况

for my $f(@files)
{
    my ($name, $type) = split/\./, $f;
    my ($genus, $species) = split(/_/, $name);
    my $suffix = "_".substr($genus, 0, 2).substr($species, 0, 3);
    # print($suffix);exit;

    if ($type eq 'gff' or $type eq 'gff3')
    {
        $type = 'gff';
    }
    else
    {
        $type = 'fasta';
    }

    open IN, $f or die;
    open OUT, ">".$f.".new" or die;
    while(<IN>)
    {
        if($type eq 'fasta')
        {
            if(/^>/)
            {
                chomp;
                s/\s.*//;
                s/$/$suffix/;
                $_.="\n";
            }
        }
        elsif($type eq 'gff')
        {
            # my $mRNA_count = `awk '\$3=="mRNA"' $f |wc -l`;
            # my $gene_count = `awk '\$3=="gene"' $f |wc -l`;
            # if($gene_count > 0 && $mRNA_count ne $gene_count)
            # {
            #     die "Unequal gene and mRNA number";
            # }
            chomp;
            if($_ eq "")
            {
                next;
            }
            my @e=split/\t/, $_;
            $e[0].=$suffix;

            if($e[2] eq 'gene')
            {
                next;
            }

            if($e[2] eq 'mRNA')
            {
                $e[-1]=~s/;Parent.*//;
            }
            my @ee = split(/;/, $e[-1]);
            for my $ee (@ee)
            {
                if($ee=~/ID/ or $ee =~/Parent/)
                {
                    $ee=~s/$/$suffix/;
                    $ee=~s/TR\.//;
                }
            }
            @e[-1] = join(';', @ee);
            $_ = join("\t", @e);
            $_ .= "\n";
        }
        print OUT $_;
    }
    close IN;
    close OUT;
}
