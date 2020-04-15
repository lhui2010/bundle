#!/usr/bin/perl -w

if(@ARGV<2)
{
        print "$0 selfblat.psl source.fa\n";
        exit;
}

open PSL, $ARGV[0] or die;
open FA, $ARGV[1] or die;

$coverage=0.3;
$minimum_match = 40;
$cluster_count = 0;
while(<PSL>)
{
        chomp;
        if($_ eq "" or /psLayout/ or /match/ or /-----/)
        {
                next;
        }


        @e=split;

        next if( $e[9] eq $e[13]);

        $match = $e[0];
        $qname = $e[9];
        $qsize = $e[10];
	$qstart = $e[11];
	$qend = $e[12];
        $tname = $e[13];
        $tsize = $e[14];
	$tstart = $e[15];
	$tend = $e[16];

#must be over hanging
	
	next unless(($tend == 0 or $tend == $tsize or $tstart == 0 or $tstart == $tsize) and ($qstart == $qsize or $qstart == 0 or $qend == $qsize or $qend == 0));

        if (($match/$qsize >= $coverage or $match/$tsize >= $coverage) and $match > $minimum_match)#overlapping threshold
        {
                my $temp = "";
                next if (exists $group{$tname} and exists $group{$qname});
                if(exists $group{$tname})
                {
                        $group{$qname} = $group{$tname};
                        push @{$cluster[$group{$qname}]}, $qname;
                }
                elsif(exists $group{$qname})
                {
                        $group{$tname} = $group{$qname};
                        push @{$cluster[$group{$tname}]}, $tname;
                }
                else
                {
                        $cluster[$cluster_count] = [$tname, $qname];
                        $group{$qname} = $cluster_count;
                        $group{$tname} = $cluster_count;
                        $cluster_count++;
                }
        }
}

$/='>';
$_=<FA>;
while(<FA>)
{
        chomp;

        (my $name, my $seq) = split /\n/, $_, 2;
        my @e=split /\s+/, $name;
        my $realname = $e[0];
        $fa{$realname} = ">$_";
}

$this_id = 0;
for $this_cluster(@cluster)
{
#       print clustered fasta
#       ..

        print "#SuperScaffold-$this_id\n";
        for (my $i = 0; $i<@{$this_cluster}; $i++)
        {
                #print $$this_cluster[$i], "\n";
                print $fa{$$this_cluster[$i]};
                delete $fa{$$this_cluster[$i]};
        }
        $this_id++;
}

#printing singletons
for $this_fa(sort keys %fa)
{
#       print singleton
#       ..
        print "#SuperScaffold-$this_id\n";
	print $fa{$this_fa};
        $this_id++;
#        print "#key\t$this_fa\n"; forget what for...
        #delete $fa{$this_fa};
}

