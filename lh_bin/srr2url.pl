#ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR066/SRR066562/SRR066562.sra

while(<>)
{
    chomp;
    my $sra=$_;
    my $pre = substr($_, 0, 6);

    print "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/$pre/$sra/$sra.sra\n";
}
