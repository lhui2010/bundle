#!/usr/bin/perl -w
#000042F|arrow_np1212	genBlastG	mRNA	1634599	1636723	141.141	-	.	ID=sp_P92994_TCMO_ARATH-R3-3-A1;Name=sp_P92994_TCMO_ARATH
#000042F|arrow_np1212	genBlastG	CDS	1635917	1636723	.	-	.	ID=sp_P92994_TCMO_ARATH-R3-3-A1-E1;Parent=sp_P92994_TCMO_ARATH-R3-3-A1
#000042F|arrow_np1212	genBlastG	CDS	1634839	1635351	.	-	.	ID=sp_P92994_TCMO_ARATH-R3-3-A1-E2;Parent=sp_P92994_TCMO_ARATH-R3-3-A1
#000042F|arrow_np1212	genBlastG	CDS	1634599	1634823	.	-	.	ID=sp_P92994_TCMO_ARATH-R3-3-A1-E3;Parent=sp_P92994_TCMO_ARATH-R3-3-A1
#000042F|arrow_np1212	genBlastG	mRNA	1153504	1156035	51.2848	-	.	ID=sp_Q9SVZ9_ORP4A_ARATH-R3-3-A1;Name=sp_Q9SVZ9_ORP4A_ARATH
#
#transcript
#coding_exon

while(<>)
{
    my @e=split;
    if($e[2] eq 'transcript')
    {
        s/transcript/mRNA/;
        my $line = $_;
        s/mRNA/gene/;
        s/;/G;/;
        print;
        $line =~ s/ID=(.*);Name=(.*)$/ID=$1;Name=$2;Parent=$1G/;
        print $line;
    }
    elsif($e[2] eq 'coding_exon')
    {
        s/coding_exon/CDS/;
        my $line = $_;
        s/CDS/exon/;
        print;
        print $line;
    }
}
