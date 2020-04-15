#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;

my ($gw_file, $pep_file, $alignRate_cutoff, $identity_cutoff, $type);
my ($Verbose,$Help);
GetOptions(
        "gw:s"=>\$gw_file,
        "pep:s"=>\$pep_file,
        "ac:s"=>\$alignRate_cutoff,
        "id:s"=>\$identity_cutoff,
        "type:s"=>\$type,
        "help"=>\$Help,
        "verbose"=>\$Verbose,
);
if (!$gw_file || !$pep_file || $Help) {
        print <<"Usage End.";

Description:
        This program is used for detect mutation(including indel) from genewise result.

        Version: 1.0  Date: 2009-11-11
        Author:  liqiye <liqiye\@genomics.org.cn>

Usage:
        --gw           genewise result
        --pep          optional, the query pep sequence for getting the length of each pep sequence. It is just used to
                       filter some result with the aligning rate when --ac is given.
        --ac           optional, the aligning rate cutoff for alignment, default=0.5
        --id           optional, the identity cut off for alignment, default=50;
        --type         numeric, 1 means output alignment info; 2 means output mutation info, default=2
        --help         output help information to screen

Example:
        perl

Usage End.

        exit;
}
$alignRate_cutoff = defined($alignRate_cutoff) ? $alignRate_cutoff : 0.5;
$identity_cutoff = defined($identity_cutoff) ? $identity_cutoff : 50;
$type = defined($type) ? $type : 2;
#################################################################################################################

## read in standard codon list for detecting site mutation.
my $codon_list = "/lustre/user/wwlab01/project/rice/bin/codon.txt";
my %codonList;
open IN, $codon_list;
while (<IN>) {
        my @info = split /\s+/;
        my $aa = shift @info;
        foreach my $codon (@info) {
                push @{$codonList{$aa}}, $codon;
        }
}
close IN;

## get pep length
my %pepLen;
&pep_parser($pep_file, \%pepLen);


## detect mutation
my @alignment;
my @result;
open IN, $gw_file;
$/ = "//\ngenewise";
while (<IN>) {
        /Query\s+protein:\s+(.+?)\n.+Target\s+Sequence\s+(.+?)\nStrand:\s+(.+?)\n/s;
        my ($query, $target, $strand) = ($1, $2, $3);
        #my ($chr, $chr_bg) = (split /_/, $target)[1,2];
#       my @chr_tmp = (split /##/, $target);
#       my ($chr, $chr_bg) = ($chr_tmp[-3], $chr_tmp[-2]);
        my $chr = $target;
        my $match = (length($query) > 15) ? substr($query, 0, 15) : $query;

        ## get alignment blocks
        my @blocks;
        if ($strand eq "both") {
                /genewise\s+output(.+)genewise\s+output/s;
                push @blocks, $1;
                /.+genewise\s+output(.+)\/\//s;
                push @blocks, $1;
        } else {
                push @blocks, $_;
        }
#print $blocks[0]."\n\n\nss\n\n".$blocks[1];
#exit;
        ## deal with each block to get mutation information
        foreach my $block (@blocks) {
		my $intron_num = 0;
		$intron_num++ while ($block =~ m/Intron/g);
#		print "$intron_num\n";
                my @info = split /\n/, $block;
                ## get alignment score
                $block =~ /Score\s+(\S+)\s+bits\s+over\s+entire\s+alignment/;
                my $score = $1;

                ## get the start site of the query pep
                my $aa_bg;
                for (my $i = 0; $i < @info; $i ++) {
                        if ($info[$i] =~ /^$match/) {
                                $aa_bg = (split /\s+/, $info[$i])[1];
                                last;
                        }
                }

                ## link the align block
                my @nuc_align_block;
                my @pep_align_block;
                my @exon_location;
                my ($query_pep, $match_pep, $target_pep, $nuc_line1, $nuc_line2, $nuc_line3);
                my $match_len;
                for (my $i = 0; $i < @info; $i ++) {
                        if ($info[$i] =~ /^$match/) {
                                ## determine line len
                                my $temp1 = $info[$i];
                                $i += 2;
                                my $temp2 = $info[$i];
                                $i += 1;
                                my $temp3 = $info[$i];
                                $i += 1;
                                my $temp4 = $info[$i];
                                $i += 1;
                                my $temp5 = $info[$i];
                                my $line_len = &determine_line_len($temp1, $temp2, $temp3, $temp4, $temp5);
                                $i -= 5;

                                ## nuc line2
                                $i += 4;
                                substr($info[$i], 0, 21) = "";
                                $info[$i] = substr($info[$i], 0, $line_len);
                                $nuc_line2 .= $info[$i];

                                ## query pep
                                $i -= 4;
                                substr($info[$i], 0, 21) = "";
                                $info[$i] = substr($info[$i], 0, $line_len);
                                $query_pep .= $info[$i];

                                ## match pep
                                $i += 1;
                                substr($info[$i], 0, 21) = "";
                                $info[$i] = substr($info[$i], 0, $line_len);
                                $match_pep .= $info[$i];

                                ## target pep
                                #$i += 2;
                                $i += 1;
                                substr($info[$i], 0, 21) = "";
                                $info[$i] = substr($info[$i], 0, $line_len);
                                $target_pep .= $info[$i];

                                ## nuc line1
                                $i += 1;
                                substr($info[$i], 0, 21) = "";
                                $info[$i] = substr($info[$i], 0, $line_len);
                                $nuc_line1 .= $info[$i];

                                ## nuc line3
                                $i += 2;
                                substr($info[$i], 0, 21) = "";
                                $info[$i] = substr($info[$i], 0, $line_len);
                                $nuc_line3 .= $info[$i];

                                $i += 2;
                        } elsif ($info[$i] =~ /\s+Exon\s+(\d+)\s+(\d+)\s+phase\s+(\d+)/) {
                                my ($exon_bg, $exon_ed, $phase) = ($1, $2, $3);
                        #       $exon_bg = $exon_bg + $chr_bg - 1;
                        #       $exon_ed = $exon_ed + $chr_bg - 1;
                                push @nuc_align_block, ($exon_bg, $exon_ed);
                                push @exon_location, [$exon_bg, $exon_ed, $phase];
                        } elsif ($info[$i] =~ /\s+Supporting\s+\d+\s+\d+\s+(\d+)\s+(\d+)/) {
                                push @pep_align_block, ($1, $2);
                                $match_len += (abs($1 - $2) + 1);
                        }
                }

                ## determine strand
                $strand = "+" if $nuc_align_block[0] < $nuc_align_block[1];
                $strand = "-" if $nuc_align_block[0] > $nuc_align_block[1];

                ## calculate identity
                $match_pep =~ s/\s+|\+//g;
                my $identity = sprintf "%.2f", length($match_pep) / $match_len * 100;
                next unless $identity > $identity_cutoff;

                ## get the start and end site of the nucleotide sequence in the alignment
                @nuc_align_block = sort {$a <=> $b} @nuc_align_block;
                my ($nuc_align_bg, $nuc_align_ed) = ($nuc_align_block[0], $nuc_align_block[-1]);

                ## get the start and end site of the query pep in the alignment
                @pep_align_block = sort {$a <=> $b} @pep_align_block;
                my ($pep_align_bg, $pep_align_ed) = ($pep_align_block[0], $pep_align_block[-1]);
                ## filter some result with too low aligning rate
                unless( defined $pepLen{$query}){
			print "###".$query."\n";
			next;
		}
		my $alignRate = sprintf "%.4f", ($pep_align_ed - $pep_align_bg + 1)/$pepLen{$query};
                next unless $alignRate > $alignRate_cutoff;
                push @alignment, "$query\t$pepLen{$query}\t$pep_align_bg\t$pep_align_ed\t$strand\t$chr\t$nuc_align_bg\t$nuc_align_ed\t$alignRate\t$score\t$identity\n";


                ## get and check the length of each line
                my $query_pep_len = length($query_pep);
                my $target_pep_len = length($target_pep);
                my $nuc_line1_len = length($nuc_line1);
                my $nuc_line2_len = length($nuc_line2);
                my $nuc_line3_len = length($nuc_line3);
                die "$query len error!\n" unless $query_pep_len == $target_pep_len && $target_pep_len == $nuc_line1_len && $nuc_line1_len == $nuc_line2_len && $nuc_line2_len == $nuc_line3_len;
                #print "$query\n$query_pep\n$target_pep\n$nuc_line1\n$nuc_line2\n$nuc_line3\n\n";

                my $exon_count = 0;
                my $bp = ($strand eq "+") ? $nuc_align_bg : $nuc_align_ed;      ## set down the the postion of each base
                my $ap = $aa_bg; ## set down the position of each animo acid of the query pep
                for (my $i = 0; $i < $query_pep_len; $i ++) {
                        my $qaa = substr($query_pep, $i, 1);
                        my $taa = substr($target_pep, $i, 1);
                        my $ba1 = substr($nuc_line1, $i, 1);
                        my $ba2 = substr($nuc_line2, $i, 1);
                        my $ba3 = substr($nuc_line3, $i, 1);

                        ## detect insertion, deletion, premature termination and frameshift.
                        my ($nuc_mut_bg, $nuc_mut_ed, $pep_mut_bg, $pep_mut_ed, $mut_info, $mut_type);
                        if ($qaa eq "-") {
                                ($nuc_mut_bg, $nuc_mut_ed) = ($ba1 =~ /\d+/) ? &get_region($strand, $bp, $ba1) : &get_region($strand, $bp, 3);
                                $mut_info = $ba1 . $ba2 . $ba3;
                                ($pep_mut_bg, $pep_mut_ed) = ($ap - 1, $ap);
                                $mut_type = "I";
                        } elsif ($taa eq "-") {
                                ($nuc_mut_bg, $nuc_mut_ed) = ($strand eq "+") ? &get_region($strand, $bp-1, 2) : &get_region($strand, $bp+1, 2);
                                $mut_info = $qaa;
                                ($pep_mut_bg, $pep_mut_ed) = ($ap, $ap);
                                $mut_type = "D";
                        } elsif ($taa eq "X" && $qaa ne "U" && $qaa ne "X") {
                                ($nuc_mut_bg, $nuc_mut_ed) = &get_region($strand, $bp, 3);
                                my $stop_codon = $ba1 . $ba2 . $ba3;
                                $mut_info = &check_mutation_site($stop_codon, \@{$codonList{$qaa}});
                                die "$query mutation info miss\n" unless $mut_info;
                                ($pep_mut_bg, $pep_mut_ed) = ($ap, $ap);
                                $mut_type = "S";
                        } elsif ($taa eq "!") {
                                die "$ba1 is not numeric\n" unless $ba1 =~ /\d+/;
                                ($nuc_mut_bg, $nuc_mut_ed) = &get_region($strand, $bp, $ba1);
                                $mut_info = "$ba1:$qaa";
                                ($pep_mut_bg, $pep_mut_ed) = ($ap, $ap);
                                $mut_type = "F";
                        }
=pod
                        elsif ($qaa =~ /\s+/ && $ba2 eq "<") {
                                my $pep_tmp = substr($target_pep, $i, 22);
                                if ($pep_tmp =~ /(\S):(\S)\[\w\w\w\]/) {
                                        my ($q_splice_aa, $t_splice_aa) = ($1, $2);
                                        if ($t_splice_aa eq "X") {
                                                my $nuc_tmp = substr($nuc_line2, $i, 24);
                                                #print "$query\t$q_splice_aa\t$t_splice_aa\t$nuc_tmp\n";
                                                die unless $nuc_tmp =~ /<(\d)-.+-\d>/;
                                                my $phase = $1;
                                                #print "$query\n" unless $phase =~ /^\d+$/;
                                                my ($lb, $rb);
                                                if ($phase == 1) {
                                                        $lb = substr($nuc_line1, $i-1, 1);
                                                        $rb = substr($nuc_line1, $i+23, 2);
                                                } elsif ($phase == 2) {
                                                        $lb = substr($nuc_line1, $i-2, 2);
                                                        $rb = substr($nuc_line1, $i+23, 1);
                                                }
                                                my $stop_codon = $lb . $rb;
                                                print "$query\t$q_splice_aa\t$t_splice_aa\t$nuc_tmp\t$stop_codon\n";
                                        }
                                }
                        }
=cut
                        #push @result, "$query\t$pepLen{$query}\t$pep_align_bg\t$pep_align_ed\t$alignRate\t$strand\t$target\t$nuc_align_bg\t$nuc_align_ed\t$score\t$identity\t$mut_type\t$pep_mut_bg\t$pep_mut_ed\t$nuc_mut_bg\t$nuc_mut_ed\t$mut_info\n" if $mut_type;
                        push @result, "$query\t$pepLen{$query}\t$pep_align_bg\t$pep_align_ed\t$alignRate\t$strand\t$chr\t$nuc_align_bg\t$nuc_align_ed\t$score\t$identity\t$mut_type\t$pep_mut_bg\t$pep_mut_ed\t$nuc_mut_bg\t$nuc_mut_ed\t$mut_info\t$intron_num\n" if $mut_type;
#			print "$intron_num\n";
                        ## deal with the query animo acid position
                        if ($qaa =~ /\s+/ && $ba2 eq "<") {
                                my $temp = substr($target_pep, $i, 22);
                                if ($temp =~ /(\S):(\S)\[\w\w\w\]/) {
                                        my ($q_splice_aa, $t_splice_aa) = ($1, $2);
                                        $ap ++ unless $q_splice_aa eq "-";
                                }
                        } elsif ($qaa =~ /\w/) {
                                $ap ++;
                        }

                        ## deal with the base position
                        if ($taa =~ /\s+/ && $ba2 eq "<") {
                                $i += 22;
                                $exon_count ++;
                                $bp = $exon_location[$exon_count]->[0];
                        } elsif ($taa eq "!") {
                                die "$ba1 is not numeric\n" unless $ba1 =~ /\d+/;
                                my $next_ba = substr($nuc_line1, $i+1, 1);
                                my $next_ba2 = substr($nuc_line2, $i+1, 1);
                                $exon_count ++ unless ($next_ba=~/[ATCGN]/ && $next_ba2!~/[atcgnN]/) || $next_ba =~ /\d+/;
                                #$bp = $exon_location[$exon_count]->[0];
                                if ($next_ba =~ /\d+/) {
                                        $bp += $ba1 if $strand eq "+";
                                        $bp -= $ba1 if $strand eq "-";
                                } else {
                                        $bp = $exon_location[$exon_count]->[0];
                                }
                        } elsif ($taa eq "-") {
                                next;
                        } else {
                                if ($ba1 =~ /[atcgn]/) {
                                        $bp += 3 if $strand eq "+";
                                        $bp -= 3 if $strand eq "-";
                                } elsif ($ba1 =~ /N/ && $ba2 =~ /[atcgnN]/ && $ba3 =~ /[atcgnN]/) {
                                        $bp += 3 if $strand eq "+";
                                        $bp -= 3 if $strand eq "-";
                                } elsif ($ba1 =~ /[ATCGN]/) {
                                        $bp += 1 if $strand eq "+";
                                        $bp -= 1 if $strand eq "-";
                                }
                        }
                }
        }
}
$/ = "\n";
close IN;

if ($type == 1) {
        #print foreach @alignment;
        my %line_to_query;
        my %line_to_score;
        foreach (@alignment) {
                my ($query, $score) = (split /\s+/)[0,-1];
                $line_to_query{$_} = $query;
                $line_to_score{$_} = $score;
        }
        print "#GeneID\toriginalLength\tpredict_bg\tpredict_ed\tstrand\tchr\tchr_bg\tchr_ed\talignRate\tscore\tidentity\n";
        foreach (sort {$line_to_query{$a} cmp $line_to_query{$b} or $line_to_score{$b} <=> $line_to_score{$a}} @alignment) {
                print;
        }

} elsif ($type == 2) {
        #print foreach @result;
        ## link colinear insertion and deletion

        for (my $i = 0; $i < @result; $i ++) {
                last unless $result[$i+1];
                my @info1 = split /\s+/, $result[$i];
                my @info2 = split /\s+/, $result[$i+1];
                next unless $info1[0] eq $info2[0];
                if ($info1[11] eq "I" && $info2[11] eq "I") {
                        my $distant = abs($info2[14] - $info1[15]);
                        if ($distant == 1) {
                                $info2[14] = $info1[14];
                                $info2[16] = $info1[16] . $info2[16];
                                $result[$i+1] = join "\t", @info2;
                                $result[$i+1] .= "\n";
                                $result[$i] = "";
                        }
                } elsif ($info1[11] eq "D" && $info2[11] eq "D") {
                        if ($info1[14] == $info2[14] && $info1[15] == $info2[15]) {
                                $info2[16] = $info1[16] . $info2[16];
                                $info2[12] = $info1[12];
                                $result[$i+1] = join "\t", @info2;
                                $result[$i+1] .= "\n";
                                $result[$i] = "";
                        }
                }
        }

        my @final_result;
        foreach (@result) {
                next unless $_;
                push @final_result, $_;
        }
        @result = ();

        print ">>query\tpepLen{query}\tpep_align_bg\tpep_align_ed\talignRate\tstrand\tchr\tnuc_align_bg\tnuc_align_ed\tscore\tmut_type\tpep_mut_bg\tpep_mut_ed\tnuc_mut_bg\tnuc_mut_ed\tmut_len\tintron_num\tmut_info\n";

my %check_final;

        foreach (@final_result) {
                my @info = split /\s+/;
                ## get the mutation length
#		print "$info[16]\t\t$info[17]\n";
                my $len;
                if ($info[11] eq "I") {
                        $len = abs($info[14] - $info[15]) + 1;
                } elsif ($info[11] eq "D") {
                        $len = length($info[16]) * 3;
                } elsif ($info[11] eq "S") {
                        my @temp = split /,/, $info[16];
                        $len = @temp;
                } elsif ($info[11] eq "F") {
                        $len = (split /:/, $info[16])[0];
                }
                ## add the mutation length to result
                push @info, $info[16];
                $info[16] = $len;
                ($info[14], $info[15]) = sort {$a <=> $b} ($info[14], $info[15]);
                #my $out = join "\t", @info;
                #print "$out\n";
#		print"$info[16]\taaaaaaaaaaaa\t$info[17]\tbbbbbbbbbbbbbbbb\t$info[18]\n";
                my $out;
                for (my $i = 0; $i < @info; $i ++) {
                        next if $i == 10;
                        $out .= "$info[$i]\t";
                }
                $out =~ s/\s+$//;

                ## correct the bug for frameshift detection.
                my @correct_result = split /\s+/, $out;
                $correct_result[10] = "F" if $correct_result[15] % 3 && $correct_result[16] !~ /=>/;
                $out = join "\t", @correct_result;
                unless ( defined $check_final{$_} ) {print "$out\n";}
		$check_final{$_}=1;

        }

}


#########################
##### subroutine ########
#########################
sub check_mutation_site {
my $query_codon = shift;
my $codons_ref = shift;
$query_codon = uc($query_codon);
my %mutation;
foreach my $ref_codon (@$codons_ref) {
        for (my $i = 0; $i < 3; $i ++) {
                my $qba = substr($query_codon, $i, 1);
                my $rba = substr($ref_codon, $i, 1);
                my $site = $i + 1;
                if ($qba ne $rba) {
                        push @{$mutation{"$ref_codon=>$query_codon"}}, "$site$rba=>$qba";
                        #print "$ref_codon\t$query_codon\t$site:$rba=>$qba\n";
                }
        }
}
my %key_to_num;
foreach my $key (keys %mutation) {
        my $num = @{$mutation{$key}};
        $key_to_num{$key} = $num;
}
foreach my $codon (sort {$key_to_num{$a} <=> $key_to_num{$b}} keys %mutation) {
        my $info = join ",", @{$mutation{$codon}};
        return "${codon}:$info";
        last;
        #print "$codon\t@{$mutation{$codon}}\n";
}
}
#########################

#########################
sub get_region {
my $strand = shift;
my $bg = shift;
my $len = shift;
my $ed;
if ($strand eq "+") {
        $ed = $bg + $len - 1;
} elsif ($strand eq "-") {
        $ed = $bg - $len + 1;
}
return ($bg, $ed);
}
#########################


#########################
sub determine_line_len {
my @lines = @_;
my @len;
foreach my $line (@lines) {
        substr($line, 0, 21) = "";
        $line =~ s/\s+$//;
        push @len, length($line);
}
@len = sort {$a <=> $b} @len;
return $len[-1];
}
#########################

#########################
sub pep_parser {
my $in_file = shift;
my $ref = shift;
open IN, $in_file;
$/ = ">";
<IN>;
$/="\n";
while (<IN>) {
	chomp;
	my $id=$1 if (/^(\S+)/);
	$/=">";
	my $line=<IN>;
	chomp$line;
	$line=~s/\*//g;
	$line=~s/\n//g;
	$line=~s/\s//g;
	$ref->{$id}=length($line);
	$/="\n";

#        /(.+)\n/;
#	print Dumper $_;
#	print Dumper "please!";
#        if(defined $1){
#	my $id = (split /\s+/, $1)[0];
#	my $id=$1 if (/parent=(\S+),/);
#        s/.+\n//;
#	s/\*//g;
#        s/\s+|>//g;
#        my $len = length($_);
#        $ref->{$id} = $len;
#	print Dumper $1;
#	print Dumper $id;
#	print Dumper $ref->{$id};
#	}
}
#$/ = "\n";
#close IN;
}
#########################



