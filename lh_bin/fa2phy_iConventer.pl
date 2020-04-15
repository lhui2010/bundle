#!/usr/bin/perl

# iConvert.pl - Convert tree files and sequence files between commun formats.
# This program supports FASTA, PIR, Clustal and Nexus sequence files and Phylip, Nexus, Clustal tree files.
#
# Written: June 4, 2007
# Author:  Mathieu Fourment, Macquarie University
#
# USAGE:   perl iconvert.pl -i inputfile -f output format [-o outputfile] [-w] [-s dna or protein]
# EXEMPLE: perl iconvert.pl -i sequences.fa -f nexus
#
# Compatible sequence file: FASTA .fa, PIR .pir, Clustal alignment .aln, Nexus .nxs, Phylip .phy
# Compatible tree file:     Phylip .phy, Clustal tree .ph, Nexus .nex
#
# DESCRIPTION
#   -i <input file>
#   -f <output format>
#      For sequence file use: fasta, nexus, clustal, pir, phylip
#      For tree file use:     nexus, clustal, phylip
#   -o <output file> or standard output if not specified
#   -s <dna|protein> Specify if dna or protein sequence when converting FROM fasta, phylip or clustal TO pir or nexus
#   -w Use \r\n for windows output

use Getopt::Long;
use warnings;
use strict;

my $numArgs = $#ARGV + 1;
my ( $help, $input, $out, $format_out, $seq_type, $windows, $branch, $delete );

GetOptions(
	'help|h'     => \$help,
	'in|i=s'     => \$input,
	'out|o=s'    => \$out,
	'format|f=s' => \$format_out,
	'type|t=s'   => \$seq_type,
	'windows|w'  => \$windows,
	'branch|b'   => \$branch,
	'delete|d'   => \$delete,
);

help() if ( $help or $numArgs == 0 );

my %input_type = (
	fasta    => "FASTA",
	pirp     => "PIR protein",
	pirn     => "PIR nucleaotide",
	clustal  => "Clustal",
	phylip   => "Phylip",
	dnexus   => "DATA nexus",
	tnexus   => "Nexus tree",
	tclustal => "Clustal tree",
	tphylip  => "Phylip tree"
);

$format_out = lc($format_out);
my $escape = "\n";
if ( defined $windows ) {
	$escape = "\r\n";
}
my @list_seq;
my $tree;

# check outfile if specified
if ( defined $out ) {
	if ( -e $out ) {
		my $a;
		do {
			print STDERR "Do you want to overwrite $out [Y/n]";
			$a = <STDIN>;
			chop($a);
		} while ( !( $a ne "y" or $a ne "n" or $a ne "" ) );
		exit(0) if ( $a eq "n" );
	}
	$out = ">" . $out;
}
else { $out = ">&STDOUT"; }

# get the format of the input file
my $format_in = test();

print STDERR "format in: $format_in\n";

if (
	(
		   $format_in eq "fasta"
		or $format_in eq "clustal"
		or $format_in eq "phylip"
	) & ( $format_out eq "nexus" or $format_out eq "pir" ) & !$seq_type
  )
{
	my $a;
	do {
		print STDERR "Type d if it is a DNA sequence file or p if it is a protein sequence file [D/p]";
		$a = <STDIN>;
		chop($a);
		$a = uc($a);
	} while ( !( $a ne "d" | $a ne "p" | $a ne "" ) );
	exit(0) if ( $a eq "n" );
	$seq_type = "DNA";
	$seq_type = "protein" if ( $a eq "p" );
}

print STDERR "\nInput format: $input_type{$format_in} file\n";
print STDERR "Output format: $format_out file\n\n";

##########################################
# READ ALIGNMENT
##########################################
readFasta()   if ( $format_in eq "fasta" );
readPir()     if ( $format_in eq "pirp" );
readPir()     if ( $format_in eq "pirn" );
readNexus()   if ( $format_in eq "dnexus" );
readClustal() if ( $format_in eq "clustal" );
readPhylip()  if ( $format_in eq "phylip" );

deleteThirdPosition() if(defined $delete);

##########################################
# READ TREE
##########################################
readNexusTree()   if ( $format_in eq "tnexus" );
readClustalTree() if ( $format_in eq "tclustal" );
readPhylipTree()  if ( $format_in eq "tphylip" );

##########################################
# WRITE ALIGNMENT
##########################################
toFasta()         if ( $format_out =~ /fasta/i );
toNexus()         if ( $format_out =~ /nexus/i & $format_in =~ /^[^t]/ );
toPir($format_in) if ( $format_out =~ /pir/i );
toClustal()       if ( $format_out =~ /clustal/i & $format_in =~ /^[^t]/ );
if ( $format_out =~ /phylip/i & $format_in =~ /^[^t]/ ) {
	if ( $format_out =~ /phyml/i ) {
		toPhylip();
	}
	elsif ( $format_out =~ /paml/i ) {
		toPhylipInterleavePAML();
	}
	elsif ( $format_out =~ /strict/i ) {
		toPhylipStrict();
	}
	else {
		toPhylip();
	}
}

##########################################
# WRITE TREE
##########################################
toNexusTree()   if ( $format_out =~ /nexus/i & $format_in   =~ /^t/ );
toPhylipTree()  if ( $format_out =~ /phylip/i & $format_in  =~ /^t/ );
toClustalTree() if ( $format_out =~ /clustal/i & $format_in =~ /^t/ );

###########################################
# DERTERMINE FILE TYPE
###########################################
sub test {
	open( IN, $input ) or die "$input $!\n";
	while (<IN>) {
		$_ =~ s/\r?\n$//g;
		last if ( $_ !~ /^$/ );
	}
	if (/^#nexus/i) {
		while (<IN>) {
			$_ =~ s/\r?\n$//g;
			return "dnexus" if (/begin data/i);
			return "tnexus" if (/begin tree/i);
		}
	}
	return "phylip"   if (/^\s*\d+\s+\d+/);
	return "tclustal" if (/^\($/);
	return "tphylip"  if (/^\(.*;$/);
	return "pirn"     if (/^>DL;/);           # PIR file nucleotide
	return "pirp"     if (/^>P1;/);           # PIR file protein
	return "fasta"    if (/^>/);              # FASTA file
	return "clustal"  if (/^CLUSTAL/i);       # FASTA file

	close(IN);
}

sub deleteThirdPosition{
	my $nchar = length( $list_seq[0]->{'seq'} );
	for(my $pos=0; $pos<length(@list_seq); $pos++) {
		my @temp = split(//,$list_seq[pos]->{'seq'});
		for(my $i=0; $i<$nchar; $i+=3){
			$temp[$i]="";
		}
		$list_seq[pos] = join("",@temp);
	}
}
###########################################
#   READ ALIGNMENT SUBROUTINES
###########################################

########
# FASTA
sub readFasta {
	my $rec = {};
	open( IN, $input ) or die "$input $!\n";
	while (<IN>) {
		s/\r?\n$//;
		next if (/^\s*$/);
		if (/^>/) {
			if ( $rec->{'seq'} ) {
				push( @list_seq, $rec );
				$rec = {};
			}
			$rec->{'name'} = substr( $_, 1 );
		}
		else {
			$rec->{'seq'} .= $_;
		}
	}
	push( @list_seq, $rec );
	close(IN);
}

########
# PHYLIP
# 10 first characters are for the sequence name, if less than 10 pad with spaces
# 60 characters per line per sequence
sub readPhylip {
	my $rec = {};
	my $pos = 0;
	my $ntax;
	my $nchar;
	open( IN, $input ) or die "$input $!\n";
	$_ = <IN>;
	if (/^\s*(\d+)\s+(\d+)/) {
		$ntax  = $1;
		$nchar = $2;
	}
	while (<IN>) {
		$_ =~ s/\r?\n$//;
		last if (/^\s*$/);
		$rec = {};
		my @a = split( /\B\s+\B/, $_ );
		if ( length(@a) > 1 ) {
			$rec->{'name'} = shift(@a);
			$rec->{'seq'} = join( "", @a );
		}
		# true phylip format
		else {
			$rec->{'seq'} = substr( $_, 10 );
			$rec->{'name'} = substr( $_, 0, 10 );
		}
		$rec->{'name'} =~ s/\s*$//;
		push( @list_seq, $rec );
	}

	$pos = 0;
	while (<IN>) {
		$_ =~ s/\r?\n$//;
		if (/^\s*$/) {
			$pos = 0;
			next;
		}
		$list_seq[$pos]->{'seq'} .= $_;
		$pos++;

	}
	close(IN);

	#check
	if ( scalar @list_seq != $ntax ) {
		print STDERR "ERROR: $ntax " . scalar @list_seq . " \n";
		exit(0);
	}
	for ( my $i = 0 ; $i < scalar @list_seq ; $i++ ) {
		$list_seq[$i]->{'seq'} =~ s/\s//;
		if ( length( $list_seq[$i]->{'seq'} ) != $nchar ) {
			print STDERR "ERROR ". $list_seq[$i]->{'name'}. " is longer or shorter than the number of characters specefied in the header: $nchar ". length( $list_seq[$i]->{'seq'} ) . "\n";
			exit(0);
		}
	}
}

########
# NEXUS
sub readNexus {
	my %s;
	my $pos = 0;
	my $rec = {};
	open( IN, $input ) or die "$input $!\n";
	$_ = <IN> while ( $_ !~ /matrix/ );
	while (<IN>) {
		$_ =~ s/\r?\n$//;
		next if (/^\s*$/);
		last if (/^\s*;/);
		my $name = substr( $_, 0, 25 );    #print STDERR $name."\n";
		if ( !$s{$name} ) {
			$pos++;
			$rec = {};
			$rec->{'seq'} = substr( $_, 27 );
			$rec->{'name'} = $name;
			push( @list_seq, $rec );
			$s{$name} = $pos;

		}
		else {
			$list_seq[ $s{$name} - 1 ]->{'seq'} .= substr( $_, 27 );
		}

	}
	close(IN);
}

########
# CLUSTAL
sub readClustal {
	my %s;
	my $pos = 0;
	my $rec = {};
	my $name;
	open( IN, $input ) or die "$input $!\n";
	<IN>;
	<IN>;
	<IN>;

	while (<IN>) {
		$_ =~ s/\r?\n$//;
		next if (/^\s*$/);
		next if (/^\s{25,}[\* ]*$/);
		last if (/^;/);
		my $name = substr( $_, 0, 25 );
		if ( !$s{$name} ) {
			$pos++;
			$rec = {};
			$rec->{'seq'} = substr( $_, 31 );
			$rec->{'name'} = $name;
			push( @list_seq, $rec );
			$s{$name} = $pos;

		}
		else {
			$list_seq[ $s{$name} - 1 ]->{'seq'} .= substr( $_, 31 );
		}

	}
	close(IN);
}

########
# PIR
sub readPir {
	my $what = @_;
	my $rec  = {};
	open( IN, $input ) or die "$input $!\n";
	while (<IN>) {
		$_ =~ s/\r?\n$//;
		next if (/^\s*$/);
		if (/^>(P1|DL);\s?(.+)/) {
			if ( $rec->{'seq'} ) {
				$rec->{'seq'} =~ s/\*?$//;
				push( @list_seq, $rec );
				$rec = {};
			}
			$rec->{'name'} = $2;
			$_ = <IN>;
		}
		else {
			$rec->{'seq'} .= $_;
		}
	}
	$rec->{'seq'} =~ s/\*?$//;
	push( @list_seq, $rec );
	close(IN);
}

#sub readNexusTree{
#  open(IN,$input) or die "$input $!\n";
#  while(<IN>){
#    $_=~s/\r?\n$//;
#    # \stree t1 = [-184.103] or ^tree PAUP_1 = [&U]
#    # anything but ( ) [ ] between the square brackets if the expression is presents
#    #$tree=$1 if(/^\s*tree \S+ = \[[^\(\)\[\]]+\] (.+)/);
#    $hash{$1}=$2 if(/^\s+(\d{1,})\s(.+[^,]),?$/);
#    $tree=$1 if(/^\s*tree[^=]*=\s*(?:\[[^\(\)\[\]]+\])?(.+)\s*/i);
#  }
#  close(IN);
#}

sub readNexusTree {
	my %taxa;
	my $nb;
	open( IN, $input ) or die "$input $!\n";
	while (<IN>) {
		$_ =~ s/\r?\n$//;

   		# \stree t1 = [-184.103] or ^tree PAUP_1 = [&U]
    	# anything but ( ) [ ] between the square brackets if the expression is presents
   		#$tree=$1 if(/^\s*tree \S+ = \[[^\(\)\[\]]+\] (.+)/);
		if ( /^\s*begin trees;/i .. /^\s*end;\s*$/ ) {
			if ( /^\s*translate/i .. /;\s*$/ ) {
				# All on one line and no numbers
				if (/^\s*translate.*;\s*$/i) {
					s/(^\s*translate\s*)|(;\s*$)//gi;
					my @names = split(/\s*,\s*/);
					$nb = 0;
					foreach my $name (@names) {
						$taxa{ ++$nb } = $name;
					}
				}
				else {
					$taxa{$1} = $2 if (/^\s+(\d{1,})\s(.+[^,]),?$/);
				}
			}

			#if(/^\s*u?tree[^=]*=\s*(?:\[[^\(\)\[\]]+\])?\*(\(.+\);)\s*$/i){
			if (/^\s*u?tree.*=\s*(?:\[[^\(\)\[\]]+\])?\s*(\(.+\);)\s*$/i) {
				$tree = $1;
				if ( scalar keys(%taxa) != 0 ) {
					$tree =~ s/([\(,])(\d+)/$1$taxa{$2}/g;
				}
				print $tree. "\n";

				#push(@trees,$tree);
			}
		}
	}
	close(IN);
}

sub readClustalTree {
	open( IN, $input ) or die "$input $!\n";
	while (<IN>) {
		$_ =~ s/\r?\n$//;
		$tree .= $_;
	}
	close(IN);
}

sub readPhylipTree {
	open( IN, $input ) or die "$input $!\n";
	while (<IN>) {
		$_ =~ s/\r?\n$//;
		$tree .= $_;
	}
	close(IN);
}

###########################################
#   WRITE
###########################################
sub toNexus {
	open( OUT, $out ) or die "$out $!\n";
	my $nchar    = length( $list_seq[0]->{'seq'} );
	my $ntax     = scalar @list_seq;
	my $datatype = "protein";
	if ( $seq_type =~ /dna/i ) { $datatype = "DNA"; }
	println("#NEXUS");
	println("BEGIN DATA;");
	println("dimensions ntax=$ntax nchar=$nchar;");
	println("format missing=?");
	println("interleave datatype=$datatype gap= -;$escape$escape");
	println("matrix");
	my $i = 0;

	while ( $i <= $nchar ) {
		foreach my $seq (@list_seq) {

			#my $name=substr($seq->{'name'},0,20);
			#my $name = substr( $seq->{'name'}, 0, 41 );
			my $name = substr( $seq->{'name'}, 0, 45 );
			$name =~ s/[,\s\(\)\|]/_/g;
			$name =~ s/_+$//;

			#$name.=" " x (23-length($name));
			#$name .= " " x ( 43 - length($name) );
			$name .= " " x ( 48 - length($name) );
			println( $name . "  " . substr( $seq->{'seq'}, $i, 50 ) );
		}
		$i += 50;
		println("") if ( $i <= $nchar );
	}
	my $space27 = " " x ( 25 + 2 );
	println( $space27 . ";" );
	println("end;");
	close(OUT);
}

sub toPhylip {
	open( OUT, $out ) or die "$out $!\n";
	my $nchar = length( $list_seq[0]->{'seq'} );
	my $ntax  = scalar @list_seq;
	println( " " . $ntax . " " . $nchar );
	my $i = 0;
	while ( $i <= $nchar ) {
		foreach my $seq (@list_seq) {
			if ( $i == 0 ) {

				#my $name=substr($seq->{'name'},0,10);
				my $name = substr( $seq->{'name'}, 0, 45 );

				#$name=substr($seq->{'name'},0,30);
				$name =~ s/[,\s]/_/g;
				$name =~ s/_+$//;

				#$name.=" " x (10-length($name));
				$name .= " " x ( 49 - length($name) );
				println( $name . substr( $seq->{'seq'}, $i, 60 ) );
			}
			else { println( substr( $seq->{'seq'}, $i, 60 ) ); }
		}
		$i += 60;
		println();
	}
	close(OUT);
}

sub toPhylipStrict {
	open( OUT, $out ) or die "$out $!\n";
	my $nchar = length( $list_seq[0]->{'seq'} );
	my $ntax  = scalar @list_seq;
	println( " " . $ntax . " " . $nchar );
	my $i = 0;
	while ( $i <= $nchar ) {
		foreach my $seq (@list_seq) {
			if ( $i == 0 ) {
				my $name = substr( $seq->{'name'}, 0, 10 );
				$name =~ s/[,\s]/_/g;
				$name .= " " x ( 10 - length($name) );
				println( $name . substr( $seq->{'seq'}, $i, 60 ) );
			}
			else { println( substr( $seq->{'seq'}, $i, 60 ) ); }
		}
		$i += 60;
		println();
	}
	close(OUT);
}

sub toPhylipInterleavePAML {
	open( OUT, $out ) or die "$out $!\n";
	my $nchar = length( $list_seq[0]->{'seq'} );
	my $ntax  = scalar @list_seq;
	println( " " . $ntax . " " . $nchar );
	my $i = 0;
	foreach my $seq (@list_seq) {
		println( $seq->{'name'} );
		while ( $i <= $nchar ) {
			println( substr( $seq->{'seq'}, $i, 60 ) );
			$i += 60;
		}
		$i = 0;
	}

=for cut
  while($i<=$nchar){
    foreach $r(@list_seq){
      if($i==0){
        #$n=substr($r->{'name'},0,10);
        #$n=substr($r->{'name'},0,49);
        $n=substr($r->{'name'},0,25);
        $n=~s/[,\s]/_/g;
        #$n.=" " x (10-length($n));
        $n.=" " x (25-length($n));
        println($n.substr($r->{'seq'},$i,60));
      }
      else{println(substr($r->{'seq'},$i,60));}
    }
    $i+=60;
    println();
  }
=cut

	close(OUT);
}

sub toFasta {
	open( OUT, $out ) or die "$out $!\n";
	my $nchar = length( $list_seq[0]->{'seq'} );
	my $i     = 0;

	foreach my $seq (@list_seq) {
		println( ">" . $seq->{'name'} );
		while ( $i <= $nchar ) {
			println( substr( $seq->{'seq'}, $i, 60 ) );
			$i += 60;
		}
		$i = 0;
	}
	close(OUT);
}

sub toPir {
	my $what = @_;
	open( OUT, $out ) or die "$out $!\n";
	my $nchar = length( $list_seq[0]->{'seq'} );
	my $i     = 0;

	foreach my $seq (@list_seq) {
		if ( $seq_type eq "dna" ) {
			println( ">DL; " . $seq->{'name'} );
			println( $seq->{'name'} );
		}
		else {
			println( ">P1; " . $seq->{'name'} );
			println( $seq->{'name'} );
		}
		while ( $i <= $nchar ) {
			if ( $i + 75 > $nchar ) {
				println( substr( $seq->{'seq'}, $i, 75 ) . "*" );
			}
			else { println( substr( $seq->{'seq'}, $i, 75 ) ); }
			$i += 75;
		}
		$i = 0;
	}
	close(OUT);
}

sub toClustal {
	my $space31 = " " x ( 25 + 6 );

	open( OUT, $out ) or die "$out $!\n";
	println("CLUSTAL W (1.82) multiple sequence alignment");
	println();
	println();
	my $nchar = length( $list_seq[0]->{'seq'} );
	my $ntax  = scalar @list_seq;
	my $i     = 0;
	while ( $i <= $nchar ) {

		foreach my $seq (@list_seq) {
			println(
				    substr( $seq->{'name'}, 0, 25 ) 
				  . "      "
				  . substr( $seq->{'seq'}, $i, 50 ) );
		}

		prints($space31);

		for ( my $j = $i ; $j < $i + 50 ; $j++ ) {
			last if ( $j >= $nchar );
			my $idd = 0;
			for ( my $k = 0 ; $k < scalar @list_seq ; $k++ ) {
				for ( my $l = $k + 1 ; $l < scalar @list_seq ; $l++ ) {
					my $a = substr( $list_seq[$k]->{'seq'}, $j, 1 );
					my $b = substr( $list_seq[$l]->{'seq'}, $j, 1 );

					#print STDERR "$a $b\n";
					$idd++ if ( $a eq $b );
				}
			}

			if   ( $idd == ( ( $ntax * $ntax ) - $ntax ) / 2 ) { prints("*"); }
			else                                               { prints(" "); }
		}
		$i += 50;
		println();
		println();
	}
	close(OUT);
}

########################################################################################
# WRITE TREES
########################################################################################
sub toNexusTree {
	my $brlens;
	my $ntax;
	my %taxa;

	$brlens = 1 if ( $tree =~ /:/ );    # branch length present
	if ( ( $tree =~ /\[/ or $tree =~ /\d:\d/ ) & defined $branch ) {
		my @s = split( /([,;\(\)])/, $tree );
		for ( my $i = 0 ; $i < scalar @s ; $i++ ) {
			if ( $s[ $i - 1 ] eq ")" and $s[$i] =~ /^(.+):(.+)$/ ) {
				$s[$i] = ":" . $2 . "[" . $1 . "]";
			}
			splice( @s, $i, 1 ) if ( $s[$i] eq "" );
		}
		$tree = join( "", @s );
	}

	my @s = split( /([,:;\(\)])/, $tree );
	for ( my $i = 0 ; $i < scalar @s ; $i++ ) {
		if ( $s[$i] eq "" ) {
			splice( @s, $i--, 1 );
		}
		else {
			$ntax++
			  if ( ( $s[ $i - 1 ] eq "(" or $s[ $i - 1 ] eq "," )
				and $s[$i] ne "(" );
		}
	}

	open( OUT, $out ) or die "$out $!\n";
	println("#NEXUS");
	println();
	println("Begin trees;");
	println("[!");
	println(">Generated from a $input_type{$format_in} file with $0");
	println(">$input");
	println("]");
	println("\tTranslate");
	my $nb = 0;

	for ( my $i = 0 ; $i < scalar @s ; $i++ ) {
		if ( $s[$i] !~ /[:,;\(\)]/ & $s[ $i - 1 ] ne ":" ) {
			$nb++;
			if   ( $nb == $ntax ) { println( "\t\t$nb " . $s[$i] ); }
			else                  { println( "\t\t$nb " . $s[$i] . "," ); }
			$taxa{ $s[$i] } = $nb;
		}
	}
	println("\t\t;");
	prints("tree t1 = ");
	for ( my $i = 0 ; $i < scalar @s ; $i++ ) {
		if (    $s[$i] ne "("
			and $s[$i] ne ")"
			and $s[$i] ne ","
			and $s[$i] ne ";"
			and $s[$i] ne ":" )
		{
			if (    $s[ $i - 1 ] ne ":"
				and $i - 1 >= 0
				and ( $s[ $i - 1 ] ne ")" and $s[ $i + 1 ] ne ":" ) )
			{
				prints( $taxa{ $s[$i] } );
			}
			else {
				prints( $s[$i] );
			}
		}
		else {
			prints( $s[$i] );
		}
	}
	println();
	println("End;");
	close(OUT);
}

sub toPhylipTree {
	if ( ( $tree =~ /\[/ or $tree =~ /\d:\d/ ) & defined $branch ) {
		my @s = split( /([,;()])/, $tree );
		for ( my $i = 0 ; $i < scalar @s ; $i++ ) {
			if ( $s[ $i - 1 ] eq ")" & $s[$i] =~ /^(.+):(.+)$/ ) {
				$s[$i] = ":" . $2 . "[" . $1 . "]";
			}
			splice( @s, $i, 1 ) if ( $s[$i] eq "" );
		}
		$tree = join( "", @s );
	}
	my @s = split( /([,:;()])/, $tree );
	for ( my $i = 0 ; $i < scalar @s ; $i++ ) {
		splice( @s, $i--, 1 ) if ( $s[$i] eq "" );
	}
	open( OUT, $out ) or die "$out $!\n";
	for ( my $i = 0 ; $i < scalar @s ; $i++ ) {
		if ( ( $s[$i] eq "(" | $s[$i] eq ")" | $s[$i] eq "," ) &
			$s[ $i + 1 ] ne ";" )
		{
			println( $s[$i] );
		}
		else { prints( $s[$i] ); }
	}
	println();
	close(OUT);
}

sub toClustalTree {
	if ( ( $tree =~ /\[/ or $tree =~ /\d:\d/ ) & defined $branch ) {
		my @s = split( /([,;()])/, $tree );
		for ( my $i = 0 ; $i < scalar @s ; $i++ ) {
			if ( $s[ $i - 1 ] eq ")" & $s[$i] =~ /^(.+):(.+)$/ ) {
				$s[$i] = ":" . $2 . "[" . $1 . "]";
			}
			splice( @s, $i, 1 ) if ( $s[$i] eq "" );
		}
		$tree = join( "", @s );
	}
	prints($tree);

	#my @s=split(/([,:;()])/,$tree);
	#for(my $i=0;$i<scalar @s;$i++){
	#  splice(@s,$i,1) if($s[$i] eq "");
	#}
	#open(OUT,$out) or die "$out $!\n";
	#for(my $i=0;$i<scalar @s;$i++){
	#  if(%hash & $s[$i+1] eq ":"){
	#    prints($hash{$s[$i]});
	#  }
	#  else {prints($s[$i]);}
	#}
	close(OUT);
	println();
}

sub println {
	my ($line) = @_;
	print OUT $line . $escape;
}

sub prints {
	my ($line) = @_;
	print OUT $line;
}

sub help {
	print STDERR "\nUSAGE:   perl iconvert.pl -i inputfile -f output format [-o outputfile] [-w] [-s dna or protein] [-b]\n";
	print STDERR "EXEMPLE: perl iconvert.pl -i sequences.fa -f nexus\n\n";
	print STDERR "Compatible sequence file: FASTA .fa, PIR .pir, Clustal alignment .aln, Nexus .nxs, Phylip .phy\n";
	print STDERR "Compatible tree file:     Phylip .phy, Clustal tree .ph, Nexus .nex\n\n";

	print STDERR "DESCRIPTION\n";
	print STDERR "  -i <input file>\n";
	print STDERR "  -o <output file> or standard output if not specified\n";
	print STDERR "  -f <output format>\n";
	print STDERR "      For sequence file use: fasta, nexus, clustal, pir, phylip, phylip-paml (interleaved), phylip-phyml(limitted to 30 letter), phylip-strict(limitted to 10 letters)\n";
	print STDERR "      For tree file use:     nexus, clustal, phylip\n";
	print STDERR "  -w Use \\r\\n for windows output\n";
	print STDERR "  -s <dna|protein> Specify if dna or protein sequence when converting FROM fasta, phylip or clustal TO pir or nexus\n\n";
	exit 1;
}

# TODO
# bootstrap values (a,b)97:0.1 or (a,b):0.1[97]
#