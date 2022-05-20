#!/usr/bin/perl -w
#usage: xx.pl genome gff 

sub formated_print
{
        $length = 70;
        my $seq_p =shift;
        if (ref($seq_p) eq 'SCALAR') {
                for ( my $pos = 0 ; $pos < length($$seq_p) ; $pos += $length ) {
                        print substr($$seq_p, $pos, $length), "\n";
                }
        }
}

#sub Complement_Reverse
sub rc
{
        my $seq_p=shift;
        if (ref($seq_p) eq 'SCALAR') { ##deal with sequence
                $$seq_p=~tr/AGCTagctMRWSYK/TCGAtcgaKYWSRM/;
                $$seq_p=reverse($$seq_p);
        }
}
#-------------------------

if(@ARGV <1)
{print "usage: xx.pl gff genome \n";exit;}

@fnames = @ARGV;

#1, load genome

open Fgff, $fnames[0] or die $!;
open Fgenome, $fnames[1] or die $!;


$_=<Fgenome>;
chomp;
$chr_no=substr $_, 1;
$chr_no=~s/\s.*//;

my %sequence;

while(<Fgenome>)
{
	chomp;
	
	if($_=~/>/)
	{
		$sequence{$chr_no} = $buffer;
		undef $buffer;
		$chr_no = substr $_, 1;
        $chr_no=~s/\s.*//;
	}
	else
	{
		$buffer.=$_;
	}
}

$sequence{$chr_no} = $buffer;
undef $buffer;

####edbug
#print ">chr01\n";
#for($i=0; $i<length($sequence{"chr01"}); $i++)
#{print substr($sequence{"chr01"}, $i, 1);
#print "\n"if(($i+1)%50==0);
#}
#exit;
#buffer is tested to be ok
###debut

my %last_loci;

while(<Fgff>)
{
	chomp;
    next if ($_ eq "" or /^#/);

#	next unless(/CDS/ or /exon/ or /cds/);
	
	my @e=split;

	unless($e[2] eq "CDS" or $e[2] eq "cds")
	{
		next;
	}

#	my @tmp=split/=/, $e[8];
	$e[8] =~ /Parent=(.*)/;
	my $gene_name = $1;
	$gene_name =~ s/;.*//;

	@tmp = split /;/, $gene_name;#incase of xx
	$gene_name = $tmp[0];
	$gene_name =~ s/[a-z]+://i;

	my $start = $e[3]-1;
	my $length = $e[4]-$start;
	my $seq_tmp = substr $sequence{$e[0]}, $start, $length;

#    print $start, "\t", $length, "\t", $e[0], "\t";
#    print $seq_tmp, "\n";exit;
	
	if($e[6] eq "-")
    {
        &rc(\$seq_tmp) ;
        if( !exists $gene_seq{$gene_name})
        {
            $gene_seq{$gene_name}.=$seq_tmp;
        }
        elsif($last_loci{$gene_name} > $start)
        {
            $gene_seq{$gene_name}.=$seq_tmp;
        }
        else
        {
            $gene_seq{$gene_name}=$seq_tmp.$gene_seq{$gene_name};
        }
    }
    else
    {
            $gene_seq{$gene_name}.=$seq_tmp;
    }
    print STDERR $gene_name, "\n";
    $last_loci{$gene_name} = $start;
}

for my $key(sort keys %gene_seq)
{
	print ">$key (".length($gene_seq{$key}).")\n";
	&formated_print (\$gene_seq{$key});
}

