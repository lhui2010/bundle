my $prev_gene="";
my $prev=0;
while(<>)
{


        unless (/exon/)
        {
                print;
                next;
        }

        /Parent=(.*?);/;
        my $this_gene = $1;
        $this_gene=~s/.*://;

        if ($this_gene ne $prev_gene)
        {
                $prev_gene = $this_gene;
                print;
                next;
        }

        my @e=split;
#SL3.0ch04       maker_ITAG      exon    2124503 2124902 .       +       .       ID=exon:Solyc04g008500.3.1.4;Parent=mRNA:Solyc04g008500.3.1
        if($prev == 0)
        {
                $prev = $e[4];
                print;
        }
        else
        {
                my $print = $_;
                $print =~s/exon/intron/;
                my $intron_start = $prev +1;
                my $intron_end = $e[3] - 1;
                $print =~s/intron\t\d+/intron\t$intron_start/;
                $print =~s/\d+\t\./$intron_start\t./;
                print $print;
                print;
        }

        $prev_gene = $this_gene;
}
