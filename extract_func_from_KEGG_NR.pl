#!/usr/bin/perl -w
use strict;

#InterproScan
#evm.model.tig00000281_pilon.61  20D84E4822E79A9A        464     HMMSmart        SM00369 Leucine-rich repeats, typical (most populate
#Blast against NR
#QueryId        QueryLen        QueryStart      QueryEnd        SubjectId       SubjectLen      SubjectStart    SubjectEnd      Identities%     Positives%      QueryCovery%    BitScore        Evalue  SubjectDescription

# a hash storing everything
my %func;

#InterproScan
open IPR, "IPRSCAN.result" or die;
while(<IPR>)
{
    chomp;
    my ($pep_id,$crc64, $pep_len, $method, $db_entry, $func_desc, $start_match, $end_match, $eval, $match_status, $date, $interpro_entry, $intropro_entry_desc, $GO) = split /\t/, $_;
    $func{$pep_id}.="IPR:$func_desc//";
}

#NR
open NR, "NR.result" or die;
while(<NR>)
{
    chomp;
    my ($QueryId, $QueryLen, $QueryStart, $QueryEnd, $SubjectId, $SubjectLen, $SubjectStart, $SubjectEnd, $Identities, $Positives, $QueryCovery, $BitScore, $Evalue, $SubjectDescription) = split /\t/, $_;
    $func{$QueryId} .= "$SubjectDescription//";
}

#KEGG
#Query_id	Query_length	Query_start	Query_end	Subject_id	Subject_length	Subject_start	Subject_end	Identity	Positive	Gap	Align_length	Score	E_value	Query_annotation	Subject_annotation
open KEGG, "KEGG.result" or die;
while(<KEGG>)
{
    chomp;
    my ($QueryId, $QueryLen, $QueryStart, $QueryEnd, $SubjectId, $SubjectLen, $SubjectStart, $SubjectEnd, $Identities, $Positives, $QueryCovery, $BitScore, $Evalue, $SubjectDescription) = split /\t/, $_;
    $func{$QueryId} .= "$SubjectDescription//";
}

for my $k(sort keys %func)
{
    print $k, "\t", $func{$k},"\n";
}

