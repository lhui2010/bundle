
#used to generate script

open LIST, $ARGV[0] or die;
#open MAIN, $ARGV[1] or die;

while(<LIST>)
{
	chomp;
	push @list, $_;
}

print "date\n";

for $spname(@list)
{
#	print "mkdir /home/SCE/liuhui/loss_gain_v2/StartCodon/$spname";

	for ($i=1; $i<10; $i++)
	{
		print "bsub -q QN_Norm perl -w    /home/SCE/liuhui/loss_gain_v2/bin/categorize_snp_for_rep.gff.pl    /home/SCE/liuhui/loss_gain_v2/soapsnp/$spname/chromosome0$i.cns.gz    /home/SCE/lvjun/IRGSP_5.0/GFF3_representative_chr/combine/RAP3_chr0$i.combine.gff3 /home/SCE/lvjun/IRGSP_5.0/ref_chr/chromosome0$i -o /home/SCE/liuhui/loss_gain_v2/RegionAll/$spname\n";
	}

        for ($i=10; $i<13; $i++)
        {
                print "bsub -q QN_Norm perl -w    /home/SCE/liuhui/loss_gain_v2/bin/categorize_snp_for_rep.gff.pl    /home/SCE/liuhui/loss_gain_v2/soapsnp/$spname/chromosome$i.cns.gz    /home/SCE/lvjun/IRGSP_5.0/GFF3_representative_chr/combine/RAP3_chr$i.combine.gff3 /home/SCE/lvjun/IRGSP_5.0/ref_chr/chromosome$i -o /home/SCE/liuhui/loss_gain_v2/RegionAll/$spname\n";
        }
}

print "date\n";
