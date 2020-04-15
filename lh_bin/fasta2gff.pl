#!/usr/bin/perl -w
 
use strict;
use Bio::SeqIO;
 
my $verbose = 0;
 
my %gene_seq;
 
## read in gff

if(@ARGV <2)
{ 
	print "$0 GFF Genome\n";
	exit;
}
 
my ($gffname, $faname) = @ARGV;
 
my %gff;
 
open GFF, $gffname or die "fail\n";
 
warn "Reading genome\n";
 
while(<GFF>){

	next if (/mRNA/);
  my ($seqid, undef, undef, $start, $end,
      undef, $strand, undef, $attrs) = split;
 
 $attrs =~ s/.*?=//;
 $attrs =~ s/;$//;
  push @{$gff{$seqid}}, [$start, $end, $attrs, $strand];
}
 
warn "OK\n";
 
 
 
## Do the fasta
 
my $seqio = Bio::SeqIO->
  new( -file => $faname,
       -format => 'fasta' )
  or die "double fail\n";
 
 
while(my $sobj = $seqio->next_seq)
{
  my $seqid = $sobj->id;
 
  unless(defined($gff{$seqid})){
    warn "no features for $seqid\n";
    next;
  }
 
  my $seq = ">".$sobj->seq;
 
my $toprint;

for(@{$gff{$seqid}})
{
    	my ($start, $end, $attrs, $strand) = @$_;
    	warn join("\t", $start, $end, $attrs), "\n"
    	if $verbose > 0;
 
 
 #   	print ">$seqid-". $attrs.
 #     	  "/$start-$end /$strand/ (". ($end-$start+1). ")\n";
 
	$toprint = substr($seq, $start, $end-$start+1);
	
	&rc (\$toprint) if($strand eq "-");
 #   &formated_print ($toprint);
 	
	$gene_seq{$seqid."-".$attrs} .= $toprint;
 	
  }
 
  #exit;
}

for my $key (sort keys %gene_seq)
{
	print ">$key(".length($gene_seq{$key}).")\n";
	&formated_print ($gene_seq{$key});
}
 
warn "OK\n";

sub formated_print
{
        my $length = 70;
        my $seq_p =shift;
                for ( my $pos = 0 ; $pos < length($seq_p) ; $pos += $length ) {
                        print substr($seq_p, $pos, $length), "\n";
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
