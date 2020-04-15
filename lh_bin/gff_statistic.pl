#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $Outdir;
GetOptions(
    "outdir:s"=>\$Outdir,
);
$Outdir ||= ".";
$Outdir =~ s/\/$//;
my $file=shift;
my %seqinfor;
my $real_gene_number;

open IN,$file || die "$!";
while(<IN>){
    if(/\tgene\t/)
    {
        $real_gene_number++;
    }
    chomp;
    my @c=split/\t/,$_;
    if(@c > 4){
        if($c[2]=~/mRNA/){
            my $id=$1 if($c[8]=~/ID=([^\s;]+)/);
            push @{$seqinfor{$id}{mRNA}},(@c);
        }elsif($c[2]=~/CDS/){
                    my $id=$1 if($c[8]=~/Parent=([^\s;]+)/);
                    push @{$seqinfor{$id}{CDS}},[@c];
            }#elsif($c[2]=~/exon/){
    #       my $id=$1 if($c[8]=~/Parent=([^\s;]+)/);
#           push @{$seqinfor{$id}{CDS}},[@c];
#       }
    }
}
close IN;

my $gene_number=0;
my $total_mRNA_length=0;
my $total_exon_number=0;
my $total_cds_length=0;
my $total_intron_length=0;

open OUT1,">$Outdir/mRNA_length.txt" || die "can't open:$!";
open (OUT2,">>$Outdir/$file.average.txt") || die "can't open:$!";
open (OUT3,">$Outdir/cds_length.txt") || die "can't open:$!";
open (OUT4,">$Outdir/exon_length.txt") || die "can't open:$!";
open (OUT5,">$Outdir/exon_number.txt") || die "can't open:$!";
open (OUT6,">$Outdir/intron_length.txt") || die "can't open:$!";

foreach my $title(keys %seqinfor){
    $gene_number++;
    my $mRNA_length=$seqinfor{$title}{mRNA}[4]-$seqinfor{$title}{mRNA}[3]+1;
    $total_mRNA_length=$total_mRNA_length+$mRNA_length;
    print OUT1 "$title\t$mRNA_length\n";
    my $cds_length=0;
    my $number=@{$seqinfor{$title}{CDS}};
    @{$seqinfor{$title}{CDS}}=sort {$a->[3]<=>$b->[3]}@{$seqinfor{$title}{CDS}};
    $total_exon_number=$total_exon_number+$number;
    for(my $i=0;$i<$number;$i++){
        my $exon_length=$seqinfor{$title}{CDS}[$i][4]-$seqinfor{$title}{CDS}[$i][3]+1;
        $cds_length=$cds_length+$exon_length;
        print OUT4 "$title\t$exon_length\n";
        if($i>=1){
            my $intron_length=$seqinfor{$title}{CDS}[$i][3]-$seqinfor{$title}{CDS}[$i-1][4]-1;
            print OUT6 "$title\t$intron_length\n";
            $total_intron_length=$total_intron_length+$intron_length;
        }
    }
    $total_cds_length=$total_cds_length+$cds_length;
    print OUT3 "$title\t$cds_length\n";
    print OUT5 "$title\t$number\n";
}
my $average_mRNA_length=$total_mRNA_length/$gene_number;
my $average_exon_length=$total_cds_length/$total_exon_number;
my $average_cds_length=$total_cds_length/$gene_number;
my $average_exon_number=$total_exon_number/$gene_number;
my $total_intron_number=$total_exon_number-$gene_number;
my $average_intron_length=$total_intron_length/$total_intron_number;
my $total_intron_length2=$total_mRNA_length-$total_cds_length;
my $average_intron_length2=$total_intron_length2/$total_intron_number;

print OUT2 "Total number of gene:\t$real_gene_number\n";
print OUT2 "Total number of mRNA:\t$gene_number\nAverage of mRNA_length:\t$average_mRNA_length\nAverage of cds_length:\t$average_cds_length\nAverage of exon_number:\t$average_exon_number\nAverage of exon_length:\t$average_exon_length\nAverage of intron_len:\t$average_intron_length\nAverage of intron_len:\t$average_intron_length2\nTotal number of exon:\t$total_exon_number\nTotal number of intron:\t$total_intron_number\nTotal intron length:\t$total_intron_length\nTotal intron length:\t$total_intron_length2\n";
        

